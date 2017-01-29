""" 
CATALOG - class structure for dealing with astronomical catalogs,
          FITS_LDAC files, and sqlite databases.

version 0.9, 2016-01-27, michael.mommert@nau.edu
"""

# Photometry Pipeline 
# Copyright (C) 2016  Michael Mommert, michael.mommert@nau.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


import os
import sys 
import numpy
import logging
import urllib.request, urllib.error, urllib.parse
import time
import sqlite3 as sql
from scipy import spatial
from astropy.io import fits
import scipy.optimize as optimization

# pipeline-related modules (only needed for testing)
import _pp_conf 

# setup logging
logging.basicConfig(filename = _pp_conf.log_filename, 
                    level    = _pp_conf.log_level,
                    format   = _pp_conf.log_formatline, 
                    datefmt  = _pp_conf.log_datefmt)


class catalog:
    def __init__(self, catalogname, display=False):
        self.data         = [] # astropy FITS_rec object
        self.catalogname  = catalogname
        self.obstime      = [None, None] # observation midtime (JD)
        self.obj          = None # header target name
        self.origin       = '' # where does the data come from?
        self.history      = '' # catalog history
        self.magsys       = '' # [AB|Vega|instrumental]
        self.fieldnames   = {} # unified field names
        self.display      = display
        
    #### data access functions

    @property
    def shape(self):
        """
        return: tuple of number of sources and fields
        """
        return (len(self.data), len(self.data.names))

    @property
    def fields(self):
        """ 
        return: array of all available fields 
        """
        return self.data.names + list(self.fieldnames.keys())


    def __getitem__(self, ident):
        """
        return: source or field
        """
        try:
            return self.data[ident]
        except KeyError:
            return self.data[self.fieldnames[ident]]


    ##### data manipulation functions

    def reject_sources_other_than(self, condition):
        """ 
        reject sources based on condition
        input: condition
        return: number of sources left
        """

        n_raw = self.shape[0]

        self.data = self.data[condition]

        logging.info('%s:reject %s sources' % 
                     (self.catalogname, n_raw-self.shape[0]))
                     
        return len(self.data)


    def reject_sources_with(self, condition):
        """ 
        reject sources based on condition
        input: condition
        return: number of sources left
        """

        n_raw = self.shape[0]
        self.data = self.data[~condition]

        logging.info('%s:reject %s sources' % 
                     (self.catalogname, n_raw-self.shape[0]))
                     
        return n_raw - len(self.data)

        
    def add_field(self, field_name, field_array, field_type='D'):
        """
        single-field wrapper for add_fields
        """
        return self.add_fields([field_name], [field_array], [field_type])


    def add_fields(self, field_names, field_arrays, field_types):
        """ 
        add fields to self.data
        input: field_names, field_arrays, field_types
        output: number of added fields
        -------------
        WARNING: self.data != self.data.columns.array
        this might lead to discrepancies in the number of sources when
        adding fields
        """

        assert len(field_names) == len(field_arrays) == len(field_types)

        new_cols = []
        for i in range(len(field_names)):
            new_cols.append(fits.Column(name=field_names[i], 
                                        format=field_types[i], 
                                        array=numpy.array(field_arrays[i])))

        # create FITS_rec from scratch, or append to existing array
        if type(self.data) == list:
            hdu = fits.BinTableHDU.from_columns(new_cols)
        else:
            old_cols = []
            for idx, field in enumerate(self.data.columns.names):
                old_cols.append(fits.Column(name=field, 
                                        format=self.data.columns[idx].format, 
                                        array=self[field]))
                # # add field with translated fieldname, if available
                # if field in self.fieldnames.values():
                #     rev_fieldnames = {val:key for key,val in 
                #                       self.fieldnames.items()}
                #     field = rev_fieldnames[field]
                #     old_cols.append(fits.Column(name=field, 
                #                         format=self.data.columns[idx].format, 
                #                         array=self[field]))

            hdu = fits.BinTableHDU.from_columns(fits.ColDefs(old_cols) +
                                                fits.ColDefs(new_cols))

        self.data = hdu.data

        return len(field_arrays)

    


    ##### data io

    ### online catalog access
        
    def download_catalog(self, ra_deg, dec_deg, rad_deg,
                         max_sources, sort=None, save_catalog=True,
                         display_progress=True):
        """ 
        download existing catalog from VIZIER server using self.catalogname
        input: ra_deg, dec_deg, rad_deg, max_sources, (display_progress), 
               (sort=['ascending', 'descending', None])
        return: number of sources downloaded
        """

        # assemble URL with target field properties to download a FITS table
        # see: http://vizier.u-strasbg.fr/doc/asu-summary.htx
        #server = 'http://vizier.cfa.harvard.edu'
        server = 'http://vizier.u-strasbg.fr'

        # use "server + "/viz-bin/asu-txt?" for text file output
        url =  server + "/viz-bin/asu-binfits?"
        url += "-source=" + _pp_conf.allcatalogs[self.catalogname]
        url += "&-c=%010.7f%+010.7f" % (ra_deg, dec_deg)
        url += "&-c.rd=%5.3f" % (rad_deg)
        # for mag uncertainties in URAT, pos uncertainties in 2MASS
        if self.catalogname.find('SDSS') == -1:
            url += "&-out.all" 
        url += "&-out.max=%d" % max_sources 
        # sort data
        if sort is not None:
            if sort == 'ascending':
                url += "&-sort=" + _pp_conf.allcatalogs_mag[self.catalogname]
            elif sort == 'descending':
                url += "&-sort=-" + _pp_conf.allcatalogs_mag[self.catalogname]
            else:
                logging.warning('not sure how to sort catalog data (%s)' % sort)
                print('not sure how to sort catalog data (%s)' % sort)

        logging.info('accessing %s on vizier: %s' % (self.catalogname, url))

        #print url
        
        ### call server
        call = urllib.request.urlopen(url)

        if not display_progress:
            download = call.read()

        else:
            download = b''
            last_status = ''
            file_size_bytes = 0
            block_size = 81920 # download block size in bytes

            while True:
                buf = call.read(block_size)

                if not buf:
                    break

                # cumulate total file size
                file_size_bytes += len(buf)
                new_status = 'downloaded %d kB' % (file_size_bytes/1024.)
                # remove last status and display new one
                print(''.join(['\r' for i in range(len(last_status))]) +
                      new_status, end='', flush=True)

                last_status = new_status
                download += buf

            print(' - done.')

        # write FITS table to file
        fitsfile = open(self.catalogname+'.fits', 'wb')
        fitsfile.write(download)
        fitsfile.close()

        # read FITS table and assign to self.data
        try:
            hdulist = fits.open(self.catalogname+'.fits', 
                                ignore_missing_end=True)
            self.data = hdulist[1].data
        except IOError:
            if self.display:
                print(self.catalogname, 'data not successfully downloaded')
            logging.error('%s data not successfully downloaded' % 
                          self.catalogname)
            return 0
        except IndexError:
            if self.display:
                print('no data available from %s' % self.catalogname)
            logging.error('no data available from %s' % self.catalogname)
            return 0

        # remove catalog file if save_catalog == False
        if not save_catalog:
            os.remove(self.catalogname+'.fits')

        self.history = '%d sources downloaded' % self.shape[0]

        # set catalog magnitude system
        self.magsystem = _pp_conf.allcatalogs_magsys[self.catalogname]

        # set unified field names
        self.fieldnames = {'ra.deg':'RAJ2000', 'dec.deg':'DEJ2000',
                           'e_ra.deg':'e_RAJ2000', 'e_dec.deg':'e_DEJ2000'}
        self.fieldnames['ident'] = {'2MASS'  : '_2MASS',
                                    'URAT-1' : 'URAT1',
                                    'SDSS-R9': 'SDSS9',
                                    'APASS9' : 'recno',
                                    'CMC15'  : 'CMC15',
                                    'PPMXL'  : 'PPMXL',
                                    'USNO-B1': 'USNO-B1.0'}[self.catalogname]

        if 'APASS' in self.catalogname:
            self.fieldnames['gmag'] = "g_mag"
            self.fieldnames['e_gmag'] = "e_g_mag"
            self.fieldnames['rmag'] = "r_mag"
            self.fieldnames['e_rmag'] = "e_r_mag"
            self.fieldnames['imag'] = "i_mag"
            self.fieldnames['e_imag'] = "e_i_mag"
            self.fieldnames['e_Vmag'] = "e_Vmag"
            self.fieldnames['e_Bmag'] = "e_Bmag"
            
        ### catalog additions and corrections

        # determine RA and DEC uncertainties for 2MASS
        if '2MASS' in self.catalogname:
            err_a = self['errMaj']
            err_b = self['errMin']
            err_theta = self['errPA']/180.*numpy.pi
            err_theta[err_theta == 0] = 1e-9
            arc_xopt = numpy.arctan(-err_b/err_a*numpy.tan(err_theta))
            ra_err  = abs(err_a*numpy.cos(arc_xopt)*numpy.cos(err_theta)- \
                          err_b*numpy.sin(arc_xopt)*numpy.sin(err_theta))
            arc_yopt = numpy.arctan(err_b/err_a* \
                                    numpy.cos(err_theta)/numpy.sin(err_theta))
            dec_err = abs(err_a*numpy.cos(arc_yopt)*numpy.sin(err_theta)+ \
                          err_b*numpy.sin(arc_yopt)*numpy.cos(err_theta))
            # ra_err, dec_err in arcsec
            self.add_field('e_RAJ2000', ra_err/3600.) 
            self.add_field('e_DEJ2000', dec_err/3600.)


        # perform correction to AB system for SDSS
        # http://www.sdss3.org/dr8/algorithms/fluxcal.php#SDSStoAB 
        if 'SDSS' in self.catalogname:
            self.data['umag'] = self.data['umag'] - 0.04
            self.data['zmag'] = self.data['zmag'] + 0.02


        logging.info('downloaded %d stars from %s' % (self.shape[0],
                                                      self.catalogname))

        hdulist.close()

        return self.shape[0]

    
    ### FITS/LDAC interface

    def read_ldac(self, filename, fits_filename=None, maxflag=None,
                  time_keyword='MIDTIMJD', exptime_keyword='EXPTIME',
                  object_keyword='OBJECT', telescope_keyword='TEL_KEYW'):
        """ 
        read in FITS_LDAC file 
        input: LDAC filename
        return: (number of sources, number of fields)
        """

        # load LDAC file
        hdulist  = fits.open(filename, ignore_missing_end=True)

        if len(hdulist) < 3:
            print(('ERROR: %s seems to be empty; check LOG file if ' + 
                   'Source Extractor ran properly') % filename)
            logging.error(('ERROR: %s seems to be empty; check LOG file if ' + 
                   'Source Extractor ran properly') % filename)
            return None
        

        # load data array
        self.data = hdulist[2].data

        # set other properties
        telescope = ''
        for line in hdulist[1].data[0][0]:
            if telescope_keyword in line:
                telescope = line.split('\'')[1]
        self.catalogname = filename
        self.origin      = '%s;%s' % (telescope.strip(), fits_filename)
        self.magsys      = 'instrumental'
        self.fieldnames  = {'ra.deg':'XWIN_WORLD', 'dec.deg':'YWIN_WORLD'}

        # reject flagged sources (if requested)
        if maxflag is not None:
            self.reject_sources_other_than(self.data['FLAGS'] <= maxflag)
            # FLAGS <= 3: allow for blending and nearby sources

        # read data from image header, if requested
        if fits_filename is not None:
            fitsheader = fits.open(fits_filename, 
                                   ignore_missing_end=True)[0].header
            self.obstime[0] = float(fitsheader[time_keyword])
            self.obstime[1] = float(fitsheader[exptime_keyword])
            self.obj        = fitsheader[object_keyword]

        logging.info('read %d sources in %d columns from LDAC file %s' % 
                     (self.shape[0], self.shape[1], filename))

        hdulist.close()

        return self.shape

    
    def write_ldac(self, filename, wcs_keys=('XWIN_WORLD', 'YWIN_WORLD'), 
                   projection_type='STG'):
        """ 
        write data in new FITS_LDAC file (mainly for use in SCAMP) 
        input: filename, ra/dec field names, projection_type
        return: number of sources written to file
        """
    
        # derive center coordinates for this field
        try:
            ctr = ((max(self[wcs_keys[0]])+min(self[wcs_keys[0]])/2.), 
                   (max(self[wcs_keys[1]])+min(self[wcs_keys[1]])/2.))
        except KeyError:
            print('ERROR: ldac_tools.write_file cannot find keys', wcs_keys)
            return None

        ### create primary header (empty)
        hdu_primary = fits.PrimaryHDU(header=fits.Header())

        pixel_scale_deg = (1,1) # fake values
        
        ### create LDAC_IMHEAD (field information)
        # header will be created automatically with binary table 
        # create header table
        data_array = \
            [("SIMPLE  =                    T / This is a FITS file\n" +
              "BITPIX  =                    0 /\n" +
              "NAXIS   =                    2 / 2D data\n" +
              "NAXIS1  =              %5d / Number of rows\n" +
              "NAXIS2  =              %5d / Number of columns\n" + 
              "EXTEND  =                    T / FITS extensions\n" +
              "EQUINOX =        2000.00000000 / Mean equinox\n" +
              "RADESYS = 'ICRS    '           / Astrometric system\n" +
              "CTYPE1  = 'RA---%s'           / WCS projection type\n" +
              "CUNIT1  = '        '           / Axis unit\n" +
              "CRVAL1  =   %15.8e / World coordinate on this axis\n" +
              "CRPIX1  =   %15.8e / Reference pixel on this axis\n" +
              "CD1_1   =   %e / Linear projection matrix\n" +
              "CD1_2   =   0.000000000000E+00 / Linear projection matrix\n" +
              "CTYPE2  = 'DEC--%s'           / WCS projection type\n" +
              "CUNIT2  = '        '           / Axis unit\n" +
              "CRVAL2  =   %15.8e / World coordinate on this axis\n" +
              "CRPIX2  =   %15.8e / Reference pixel on this axis\n" +
              "CD2_1   =   0.000000000000E+00 / Linear projection matrix\n" +
              "CD2_2   =   %e / Linear projection matrix\n") %
             (self.shape[0], self.shape[0], \
              projection_type, ctr[0], self.shape[0]/2., pixel_scale_deg[0], \
              projection_type, ctr[1], self.shape[0]/2., pixel_scale_deg[1])]
        ### todo: check if all these header information is necessary
        ### also check if the fake pixel resolution works with scamp 

        column = fits.Column(name='Field Header Card', format='1680A', 
                             array=data_array)
        cols = fits.ColDefs([column])
        hdu_imhead = fits.BinTableHDU.from_columns(cols)
        hdu_imhead.header['EXTNAME'] = ('LDAC_IMHEAD') 
        # manually add table name to the header 

        ### create LDAC_OBJECTS (catalog data)
        # based on Source Extractor field names
        column_array = []
        nsrc = 0
        for key in self.fields:
            # world coordinates and other angles
            if key.find('WORLD') > -1 or key.find('THETA'):
                column_array.append(fits.Column \
                                    (name=key, format='1D', \
                                     unit='deg', disp='E15', \
                                     array=self[key]))
            # image coordinates and angles
            elif key.find('IMAGE') > -1:
                column_array.append(fits.Column \
                                    (name=key, format='1D', \
                                     unit='pix', disp='E12', \
                                     array=self[key]))
            # fluxes
            elif key.find('FLUX') > -1:
                column_array.append(fits.Column \
                                    (name=key, format='1D', \
                                     unit='flux', disp='E12', \
                                     array=self[key]))
            # magnitudes
            elif key.find('MAG') > -1:
                column_array.append(fits.Column \
                                    (name=key, format='1D', \
                                     unit='mag', disp='E12', \
                                     array=self[key]))
            else:
                logging.warning(('cannot interpret field %s in writing' +\
                                 'catalog %s') % (key, self.catalogname))

            nsrc = len(self[key])
                
        columns = fits.ColDefs(column_array)
        hdu_objects = fits.BinTableHDU.from_columns(columns)
        hdu_objects.header['EXTNAME'] = ('LDAC_OBJECTS')
        # manually add table name to the header 

        # combine HDUs and write file
        hdulist = fits.HDUList([hdu_primary, hdu_imhead, hdu_objects])
        hdulist.writeto(filename, clobber=True)

        logging.info('wrote %d sources from %s to LDAC file' % 
                     (nsrc, self.filename))

        return nsrc


    ### ascii interface

    def write_ascii(self, filename):
        """ 
        write catalog to ascii table
        input: target filename
        return: number of sources written to file
        """

        # prepare headerline and formatline
        legend, headerline, formatline  = '', '', ''
        for idx in range(len(self.data.names)):
            legend += '%d - %s\n' % (idx, self.data.names[idx])
            headerline += ('|%13s ' % str(self.data.names[idx]))
            if 'E' in str(self.data.formats[idx]):
                formatline += '%15E'
            elif 'D' in str(self.data.formats[idx]):
                formatline += '%15f'
            elif 'I' in str(self.data.formats[idx]) or \
                 'J' in str(self.data.formats[idx]):
                formatline += '%15d'

        # write data into file
        numpy.savetxt(filename, self.data, fmt=formatline, 
                      header=legend+headerline)

        logging.info('wrote %d sources from %s to ASCII file %s' % 
                     (self.shape[0], self.catalogname, filename))
        
        return self.shape[0]
    
    
    ### SQLite interface

    def write_database (self, filename):
        """ 
        write catalog object to SQLite database file
        input: target filename
        output: number of sources written to file
        """

        # open database file (delete existing ones)
        os.remove(filename) if os.path.exists(filename) else None
        db_conn = sql.connect(filename)
        db = db_conn.cursor()

        # create header table
        db.execute("CREATE TABLE header (" + \
                   "[name] TEXT, [origin] TEXT, [description] TEXT, " + \
                   "[magsys] TEXT, [obstime] REAL, [exptime] REAL, [obj] TEXT)")
        db.execute("INSERT INTO header VALUES (?,?,?,?,?,?,?)", 
                   (self.catalogname, self.origin, self.history,
                    self.magsys, self.obstime[0], self.obstime[1], 
                    self.obj))


        # create data table
        table_cmd = "CREATE TABLE data ("
        for key_idx, key in enumerate(self.data.names):
            db_key = key
            # translate unified fieldnames
            for unikey, unival in self.fieldnames.items():
                if key == unival:
                    db_key = unikey
            if type(self.data[key][0]) == numpy.float32 \
               or type(self.data[key][0]) == numpy.float64:
                table_cmd += "'%s' REAL" % db_key
            if type(self.data[key][0]) == numpy.int16 \
               or type(self.data[key][0]) == numpy.int32:
                table_cmd += "'%s' INTEGER" % db_key 
            if type(self.data[key][0]) == numpy.string_:
                table_cmd += "'%s' TEXT" % db_key 
            if key_idx < len(self.data.names)-1:
                table_cmd += ", "

        table_cmd += ")"
        db.execute(table_cmd)

        # create a data array in which data types are converted to SQL types
        sqltypes = {numpy.float32:numpy.float64, numpy.float64:numpy.float64,
                    numpy.int16:numpy.int64, numpy.int32:numpy.int64}
        data_cols = [self.data[key].astype(sqltypes[type(self.data[key][0])]) \
                     for key in self.data.names]
        data = [[data_cols[j][i] for j in range(len(data_cols))] \
                for i in range(len(data_cols[0]))]
        db.executemany("INSERT INTO data VALUES (" + \
                       ','.join(['?' for i in range(len(self.data.names))]) + \
                       ')', data)
        
        db_conn.commit()

        # return number of objects written to database
        db.execute("SELECT COUNT(DISTINCT %s) FROM data" % self.data.names[0])
        n_obj = db.fetchall()[0][0]
        
        logging.info('wrote %d sources from catalog %s to database file %s' %
                     (n_obj, " | ".join([self.catalogname, self.origin, 
                                         self.history]),
                      filename))

        return n_obj

    
    def read_database (self, filename):
        """ read in photometry database into catalog """

        # open database file
        try:
            db_conn = sql.connect(filename)
            db = db_conn.cursor()
        except:
            if self.display:
                print('ERROR: could not find database', filename)
                logging.error('ERROR: could not find database', filename)
            return []

        # query database header
        db.execute("SELECT * FROM header")
        rows = db.fetchall()

        self.catalogname = rows[0][0].encode('ascii')
        self.origin      = rows[0][1].encode('ascii')
        self.history     = rows[0][2].encode('ascii')
        self.magsys      = rows[0][3].encode('ascii')
        self.obstime[0]  = rows[0][4]
        self.obstime[1]  = rows[0][5]
        self.obj         = rows[0][6].encode('ascii')
        self.fieldnames  = {}

        # query database sources
        db.execute("SELECT * FROM data")
        rows = db.fetchall()
        
        # read in field names and types
        fieldnames, types = [], []
        type_dict = {float:numpy.float64, int:numpy.int64, str:numpy.string_}
        for key_idx, key in enumerate(db.description):
            fieldnames.append(key[0])
            types.append(type_dict[type(rows[0][key_idx])])

        # read in data in FITS_rec structure by creating a
        # temporary FITS table
        data = [[] for i in range(len(fieldnames))]
        for row in rows:
            for col in range(len(row)):
                data[col].append(types[col](row[col]))
        type_dict = {numpy.float64:'E', numpy.int64:'I', numpy.string_:'A'}
        columns = [fits.Column(name=fieldnames[idx],
                               format=type_dict[types[idx]],
                               array=data[idx])
                   for idx in range(len(data))]
        tbhdu = fits.BinTableHDU.from_columns(columns)

        self.data = tbhdu.data

        return self.shape[0]
            
       

    ##### filter transformations

    def lin_func(self, x, a, b):
        return a*x + b
    
    def transform_filters (self, targetfilter):
        """
        transform a given catalog into a different filter band; crop the
        resulting catalog to only those sources that have transformed magnitudes
        transformed magnitudes start with an asterisk
        input: targetfilter name
        return: number of transformed magnitudes
        """

        if len(self.data) == 0:
            return 0

        ### SDSS to BVRI
        if ('SDSS' in self.catalogname) and \
           (targetfilter in {'B', 'V', 'R', 'I'}) and \
           (self.magsystem == 'AB'):

            logging.info(('trying to transform %d SDSS sources to ' \
                          + '%s') % (self.shape[0], targetfilter))

            ### transformations based on Chonis & Gaskell 2008, AJ, 135 
            mags = numpy.array([self['gmag'], self['rmag'],
                                self['imag'],
                                self['e_gmag'],
                                self['e_rmag'],
                                self['e_imag'],
                                self['umag'],
                                self['e_umag']])

            # lbl   = {'_Bmag':0 , '_e_Bmag': 1, '_Vmag': 2, '_e_Vmag': 3, 
            #          '_Rmag': 4, '_e_Rmag': 5, '_Imag': 6, '_e_Imag': 7}
            # nmags = [numpy.zeros(self.shape[0]) for i in range(len(lbl))]

            # sort out sources that do not meet the C&G requirements
            keep_idc = (mags[1]-mags[2] > 0.08) & (mags[1]-mags[2] < 0.5) & \
                       (mags[0]-mags[1] > 0.2)  & (mags[0]-mags[1] < 1.4) & \
                       (mags[0] >= 14.5)        & (mags[0] < 19.5) & \
                       (mags[1] >= 14.5)        & (mags[1] < 19.5) & \
                       (mags[2] >= 14.5)        & (mags[2] < 19.5)
            filtered_mags = numpy.array([mags[i][keep_idc] 
                                         for i in range(len(mags))])

            # ... derive a linear best fit and remove outliers (>3 sigma)
            ri = numpy.array(filtered_mags[1]) - numpy.array(filtered_mags[2])
            gr = numpy.array(filtered_mags[0]) - numpy.array(filtered_mags[1])

            if len(ri) == 0 or len(gr) == 0:
                logging.warning('no suitable stars for transformation to %s' %
                                targetfilter)
                return 0

            param = optimization.curve_fit(self.lin_func, ri, gr, [1,0])[0]
            resid = numpy.sqrt(((ri+param[0]*gr-param[0]*param[1])/
                                (param[0]**2+1))**2+
                               (param[0]*(ri+param[0]*gr-param[0]*param[1])/
                                (param[0]**2+1)+param[1]-gr)**2)
            remove = numpy.where(numpy.array(resid) > 3.*numpy.std(resid))[0]

            keep = [idx for idx in numpy.arange(len(keep_idc))[keep_idc]
                    if idx not in set(remove)]

            # transformed magnitudes; uncertainties through Gaussian and C&G2008
            nmags = numpy.array([numpy.empty(len(mags[0])), 
                                 numpy.empty(len(mags[0]))])
            if targetfilter == 'U':
                nmags[0] = mags[6] - 0.854
                nmags[1] = numpy.sqrt(mags[7]**2 + 0.007**2)
            elif targetfilter == 'B':
                nmags[0]   = mags[0] + 0.327*(mags[0] - mags[1]) + 0.216
                nmags[1] = numpy.sqrt(((1+0.327)*mags[3])**2 \
                                      + (0.327*mags[4])**2 \
                                      + ((mags[0]-mags[1])*0.047)**2\
                                      + 0.027**2)
            elif targetfilter == 'V':
                nmags[0]   = mags[0] - 0.587*(mags[0] - mags[1]) - 0.011 
                nmags[1] = numpy.sqrt(((1+0.587)*mags[3])**2 \
                                      + (0.587*mags[4])**2 \
                                      + ((mags[0]-mags[1])*0.022)**2 \
                                      + 0.011**2)
            elif targetfilter == 'R':
                nmags[0]   = mags[1] - 0.272*(mags[1] - mags[2]) - 0.159 
                nmags[1] = numpy.sqrt(((1-0.272)*mags[4])**2 \
                                      + (0.272*mags[5])**2 \
                                      + ((mags[1]-mags[2])*0.092)**2 \
                                      + 0.022**2)
            # elif targetfilter == 'I':                    
            #     nmags[0]   = mags[2] - 0.337*(mags[1] - mags[2]) - 0.370 
            #     nmags[1] = numpy.sqrt(((1+0.337)*mags[5])**2 \
            #                           + (0.337*mags[4])**2 \
            #                           + ((mags[1]-mags[2])*0.191)**2 \
            #                           + 0.041**2)
            
            ##### experimental!!!!! Based on Connor's transformations
            elif targetfilter == 'I':                    
                nmags[0]   = mags[2] - 0.2526*(mags[1] - mags[2]) - 0.3636 
                nmags[1] = numpy.sqrt(((1+0.2526)*mags[5])**2 \
                                      + (0.2526*mags[4])**2 \
                                      + ((mags[1]-mags[2])*0.191)**2 \
                                      + 0.028**2)

            # add new filter and according uncertainty to catalog
            self.add_field('_'+targetfilter+'mag', nmags[0])
            self.add_field('_e_'+targetfilter+'mag', nmags[1])

            # get rid of sources that have not been transformed
            self.data = self.data[keep]            

            self.catalogname  += '_transformed'
            self.history += ', %d transformed to %s (Vega)' % \
                            (self.shape[0], targetfilter)
            self.magsystem = 'AB (ugriz), Vega (%s)' % targetfilter

            logging.info(('%d sources sucessfully transformed to %s' % 
                          (self.shape[0], targetfilter)))


            return self.shape[0]


        ### URAT/APASS to R/I
        # todo: include g magnitudes and account for color
        elif ('URAT' in self.catalogname or
              'APASS' in self.catalogname and 
              (targetfilter == 'R' or targetfilter == 'I')
              and self.magsystem == 'Vega'):

            logging.info('trying to transform %d %s sources to %s' % 
                          (self.shape[0], self.catalogname, targetfilter))


            ### transformations based on Chonis & Gaskell 2008, AJ, 135 
            mags = numpy.array([self['rmag'],
                                self['imag'],
                                self['e_rmag'],
                                self['e_imag']])

            # sort out sources that do not meet the C&G requirements
            keep_idc = (mags[0]-mags[1] > 0.08) & (mags[0]-mags[1] < 0.5)

            # transformed magnitudes; uncertainties through Gaussian and C&G2008
            nmags = numpy.array([numpy.empty(len(mags[0])), 
                                 numpy.empty(len(mags[0]))])

            if targetfilter == 'R':
                nmags[0]   = mags[0] - 0.272*(mags[0] - mags[1]) - 0.159 
                nmags[1] = numpy.sqrt(((1-0.272)*mags[2])**2 \
                                      + (0.272*mags[3])**2 \
                                      + ((mags[0]-mags[1])*0.092)**2 \
                                      + 0.022**2)
            elif targetfilter == 'I':                    
                nmags[0]   = mags[1] - 0.337*(mags[0] - mags[1]) - 0.370 
                nmags[1] = numpy.sqrt(((1+0.337)*mags[3])**2 \
                                      + (0.337*mags[2])**2 \
                                      + ((mags[0]-mags[1])*0.191)**2 \
                                      + 0.041**2)

            # add new filter and according uncertainty to catalog
            self.add_field('_'+targetfilter+'mag', nmags[0])
            self.add_field('_e_'+targetfilter+'mag', nmags[1])

            # get rid of sources that have not been transformed
            self.data = self.data[keep_idc]            

            self.catalogname  += '_transformed'
            self.history += ', %d transformed to %s (Vega)' % \
                            (self.shape[0], targetfilter)
            self.magsystem = 'Vega'
            
            logging.info('%d sources sucessfully transformed to %s' % \
                         (self.shape[0], targetfilter))
            
            return self.shape[0]

       
        ### 2MASS to Warner BVRI (not accounting for galactic extinction)
        elif (self.catalogname == '2MASS') and \
           (targetfilter in {'B', 'V', 'R', 'I'}) and \
           (self.magsystem == 'Vega'):

            logging.info(('trying to transform %d 2MASS sources to ' \
                          + 'Warner %s') % (self.shape[0], targetfilter))

            # transformations using the recipe by Warner 2007, MPBu
            mags  = [self['Jmag'], self['Kmag'], self['e_Jmag']]
            lbl   = {'_Bmag':0 , '_e_Bmag': 1, '_Vmag': 2, '_e_Vmag': 3, 
                     '_Rmag': 4, '_e_Rmag': 5, '_Imag': 6, '_e_Imag': 7}
            nmags = [numpy.zeros(self.shape[0]) for i in range(len(lbl))]

            for idx in range(self.shape[0]):
                keep = True
                # remove objects with extreme J-Ks color index
                cidx = mags[0][idx]-mags[1][idx]
                if cidx < -0.1 or cidx > 1.0:
                    keep = False
                # reject faint stars based on Figure 6 by Hodgkin et
                # al. 2009, MNRAS
                if mags[0][idx] > 18 or mags[1][idx] > 17:
                    keep = False

                if keep:
                    nmags[lbl['_Bmag']][idx]   = mags[0][idx] + 1.7495*cidx**3 \
                                                - 2.7785*cidx**2 \
                                                + 5.215*cidx + 0.1980
                    nmags[lbl['_e_Bmag']][idx] = numpy.sqrt(0.08*0.08 +
                                                           mags[2][idx]**2)

                    nmags[lbl['_Vmag']][idx]   = mags[0][idx] + 1.4688*cidx**3 \
                                                - 2.3250*cidx**2 \
                                                + 3.5143*cidx + 0.1496
                    nmags[lbl['_e_Vmag']][idx] = numpy.sqrt(0.05*0.05 + 
                                                           mags[2][idx]**2)

                    nmags[lbl['_Rmag']][idx]   = mags[0][idx] + 1.1230*cidx**3 \
                                                - 1.7849*cidx**2 \
                                                + 2.5105*cidx + 0.1045
                    nmags[lbl['_e_Rmag']][idx] = numpy.sqrt(0.08*0.08 + 
                                                           mags[2][idx]**2)
                
                    nmags[lbl['_Imag']][idx]   = mags[0][idx] + 0.2963*cidx**3 \
                                                - 0.4866*cidx**2 \
                                                + 1.2816*cidx + 0.0724
                    nmags[lbl['_e_Imag']][idx] = numpy.sqrt(0.03*0.03 + 
                                                           mags[2][idx]**2)
                else:
                    for mag in lbl.keys():
                        nmags[lbl[mag]][idx] = 99
              
            # append nmags arrays to catalog
            for key, idx in lbl.items():
                self.add_field(key, nmags[idx])
                
            # get rid of sources that have not been transformed
            self.data = self.data[self['_Vmag'] < 99]

            self.catalogname  += '_transformed'
            self.history += ', %d transformed to Warner %s (Vega)' % \
                            (self.shape[0], targetfilter)
            self.magsystem = 'Vega'

            logging.info('%d sources sucessfully transformed to Warner %s' % \
                         (self.shape[0], targetfilter))
            
            return self.shape[0]


        ### 2MASS to UKIRT YZJHK)
        elif (self.catalogname.find('2MASS') > -1) and \
           (targetfilter in {'Y_UKIRT', 'Z_UKIRT', 'J_UKIRT', \
                             'H_UKIRT', 'K_UKIRT'}) and \
           (self.magsystem == 'Vega'):

            logging.info(('trying to transform %d 2MASS sources to ' \
                          +'UKIRT YZJHK') % self.shape[0])

            # transformations using the recipe by Hodgkin et al. 2009, MNRAS
            mags  = [self['Jmag'], self['Hmag'], self['Kmag'],
                     self['e_Jmag'], self['e_Hmag'], 
                     self['e_Kmag']]
            lbl   = {'_Ymag': 0, '_e_Ymag': 1, '_Zmag': 2, '_e_Zmag': 3, 
                     '_Jmag': 4, '_e_Jmag': 5, '_Hmag': 6, '_e_Hmag': 7, 
                     '_Kmag': 8, '_e_Kmag': 9}
            nmags = [numpy.zeros(self.shape[0]) for i in range(len(lbl))]


            for idx in range(self.shape[0]):
                keep = True
                # remove objects with extreme J-Ks color index
                cidx = mags[0][idx]-mags[1][idx]
                if cidx < -0.1 or cidx > 1.0:
                    keep = False
                # reject faint stars based on Figure 6 by Hodgkin et
                # al. 2009, MNRAS
                if mags[0][idx] > 18 or mags[1][idx] > 17:
                    keep = False

                if keep:
                    # 0.064 (sig: 0.035) is a systematic offset
                    # between UKIRT_Z and SDSS_Z
                    nmags[lbl['_Zmag']][idx]   = mags[0][idx] \
                                                 + 0.95*(mags[0][idx]
                                                         - mags[1][idx]) + 0.064
                    nmags[lbl['_e_Zmag']][idx] = numpy.sqrt(mags[3][idx]**2 \
                                                            + 0.035**2)

                    nmags[lbl['_Ymag']][idx]   = mags[0][idx] \
                                                 + 0.5*(mags[0][idx] \
                                                        - mags[1][idx]) + 0.08
                    nmags[lbl['_e_Ymag']][idx] = mags[3][idx]

                    oldH = mags[1][idx]
                    nmags[lbl['_Hmag']][idx]   = mags[1][idx] \
                                                 + 0.07*(mags[0][idx] - \
                                                         mags[2][idx]) - 0.03
                    nmags[lbl['_e_Hmag']][idx] = mags[4][idx]

                    nmags[lbl['_Jmag']][idx]   = mags[0][idx] \
                                                 - 0.065*(mags[0][idx] - oldH)
                    nmags[lbl['_e_Jmag']][idx] = mags[3][idx]
                
                    nmags[lbl['_Kmag']][idx]    = mags[2][idx] \
                                                  + 0.01*(mags[0][idx] \
                                                          - mags[2][idx])
                    nmags[lbl['_e_Kmag']][idx]  = mags[5][idx]

                else:
                    for mag in lbl.keys():
                        nmags[lbl[mag]][idx] = 99
              
            # append nmags arrays to catalog
            for key, idx in lbl.items():
                self.add_field(key, nmags[idx])
                
            # get rid of sources that have not been transformed
            self.data = self.data[self['_Zmag'] < 99]

            self.catalogname  += '_transformed'
            self.history += ', %d transformed to UKIRT YZJHK (Vega)' %\
                            self.shape[0]
            self.magsystem = 'Vega'

            logging.info('%d sources sucessfully transformed to UKIRT YZJHK' % \
                         self.shape[0])
            
            return self.shape[0]
        

        ### SDSS to UKIRT Z
        elif (self.catalogname.find('SDSS') > -1) and \
             targetfilter == 'Z_UKIRT' and \
             self.catalogname == 'AB':

            logging.info(('trying to transform %d SDSS sources to ' \
                          + 'UKIRT Z') % self.shape[0])

            ### transformations using the recipe by Hodgkin et al. 2009, MNRAS
            # select sources for input catalog
            mags  = [self['zmag'], self['imag'], self['e_zmag']]
            lbl   = {'_Zmag': 0, '_e_Zmag': 1}
            nmags = [numpy.zeros(self.shape[0]) for i in range(len(lbl))]

            for idx in range(self.shape[0]):
                # transformations according to Hewett et al. 2006, MNRAS
                # (slope is average of III and V classes)
                nmags[lbl['_Zmag']][idx]   = mags[0][idx] - 0.01 + 0.06* \
                                             (mags[1][idx]-mags[0][idx])-0.528
                nmags[lbl['_e_Zmag']][idx] = mags[2][idx]
              
            # append nmags arrays to catalog
            for key, idx in lbl.items():
                self.add_field(key, nmags[idx])
                
            # get rid of sources that have not been transformed
            self.data = self.data[self['_Zmag'] < 99]

            self.catalogname  += '_transformed'
            self.history += ', %d transformed to UKIRT Z (Vega)' % \
                            self.shape[0]
            self.magsystem = 'AB (ugriz), Z_UKIRT (Vega)'

            logging.info('%d sources sucessfully transformed to UKIRT Z' % \
                         self.shape[0])
            
            return self.shape[0]

        else:
            if self.display:
                print('ERROR: no transformation from %s to %s available' % \
                    (self.catalogname, targetfilter))
            return 0
        

    ##### catalog operations


    def match_with(self, catalog,
                   match_keys_this_catalog=['ra.deg', 'dec.deg'],
                   match_keys_other_catalog=['ra.deg', 'dec.deg'],
                   extract_this_catalog=['ra.deg', 'dec.deg'],
                   extract_other_catalog=['ra.deg', 'dec.deg'],
                   tolerance=0.5/3600.):
        """ match sources from different catalogs based on two fields 
            (e.g., postions)
            return: requested fields for matched sources 
            note: will only match exclusive pairs
        """

        # replace unified fieldnames with actual field names
        # where applicable
        for lst in [match_keys_this_catalog, extract_this_catalog]:
            for idx, item in enumerate(lst):
                if item in self.fieldnames:
                    lst[idx] = self.fieldnames[item]
        for lst in [match_keys_other_catalog, extract_other_catalog]:
            for idx, item in enumerate(lst):
                if item in self.fieldnames:
                    lst[idx] = catalog.fieldnames[item]

        this_tree = spatial.KDTree(list(zip(self[match_keys_this_catalog[0]],
                                   self[match_keys_this_catalog[1]])))

        # kd-tree matching
        if tolerance is not None:
            other_tree = spatial.KDTree(\
                            list(zip(catalog[match_keys_other_catalog[0]],
                                catalog[match_keys_other_catalog[1]])))
            match = this_tree.query_ball_tree(other_tree, tolerance)

            indices_this_catalog  = [x for x in range(len(match)) if len(match[x]) == 1]
            indices_other_catalog = [x[0] for x in filter((lambda x: len(x)==1), match)]
            
        else:
            # will find the closest match for each target in this catalog
            other_cat = list(zip(catalog[match_keys_other_catalog[0]],
                                 catalog[match_keys_other_catalog[1]]))

            match = this_tree.query(other_cat) 

            indices_this_catalog, indices_other_catalog = [], []

            for target_idx in range(self.shape[0]):
                other_cat_indices = numpy.where(match[1]==target_idx)[0]
                if len(other_cat_indices) > 0:
                    # find closest match
                    min_idx = other_cat_indices[
                        numpy.argmin([match[0][i] for i in other_cat_indices])
                    ]

                    indices_this_catalog.append(target_idx)
                    indices_other_catalog.append(min_idx)


        ### match outputs based on indices provided and require
        ### extract_fields to be filled
        assert len(indices_this_catalog) == len(indices_other_catalog)
        indices = list(zip(indices_this_catalog, indices_other_catalog))
        # check if element is either not nan or not a float
        check_not_nan = lambda x: not numpy.isnan(x) if \
                                  (type(x) is numpy.float_) else True

        indices = filter(lambda i: all([check_not_nan(self[i[0]][key]) 
                                        for key in extract_this_catalog] \
                                       + [check_not_nan(catalog[i[1]][key]) 
                                          for key in extract_other_catalog]), 
                         indices)
        output_this_catalog  = [self[key][[i[0] for i in indices]] 
                                for key in extract_this_catalog]
        output_other_catalog = [catalog[key][[i[1] for i in indices]] 
                                for key in extract_other_catalog]

        return [output_this_catalog, output_other_catalog]





##### Test Routines

### LDAC routines

# test = catalog('ldac test')
# read = test.read_fits('testdata/test.ldac', maxflag=0)
# print '%d sources and %d columns read from ldac file' % \
#     (read[0], read[1])

# print test.write_ascii('testdata/test.ldac.dat'), \
#     'sources sucessfully written to ascii file'

# print test.write_fits('testdata/test_written.ldac'), \
#     'sources successfully written to fits file'

# print test.write_database('testdata/test.ldac.db')


### catalog download and manipulation

# cat1 = catalog('URAT-1')
# print cat1.download_catalog(80, 0, 0.5, 10000), 'sources grabbed from', cat1.catalogname
# print cat1.fields
                                                   
# cat2 = catalog('2MASS')
# print cat2.download_catalog(80, 0, 0.1, 10000), 'sources grabbed from', cat2.catalogname
# print cat2[305]
# print cat2.fields

# cat3 = catalog('SDSS-R9')
# print cat3.download_catalog(80, 40, 0.05, 10000, sort='ascending'), 'sources grabbed from', cat3.catalogname
# #print cat3[305]
# #print cat3.fields

# cat4 = catalog('USNO-B1')
# print cat4.download_catalog(80, 0, 0.5, 10000), 'sources grabbed from', cat4.catalogname
# print cat4[305]
# print cat4.fields


# print cat2.shape, cat2.history
# print cat2.transform_filters('V'), 'sources transformed to V'      2MASS to BVRI
# print cat2.shape, cat2.history

# print cat2.shape, cat2.history
# print cat2.transform_filters('K_UKIRT'), 'sources transformed to K'      #2MASS to UKIRT
# print cat2.shape, cat2.history

# print cat3.shape, cat3.history
# print cat3.transform_filters('Z_UKIRT'), 'sources transformed to Z'      #SDSS to UKIRT
# print cat3.shape, cat3.history

# print cat1.shape, cat1.history
# print cat1.transform_filters('I'), 'sources transformed to I'      #SDSS to BVRI
# print cat1.shape, cat1.history
# for i in cat1:
#    print i['_Imag'], i['rmag'], i['imag']



# print test.write_database('test.db'), 'sources written to database file'

# cat4 = catalog('')
# print cat4.read_database('test.db'), 'sources read from database file'
