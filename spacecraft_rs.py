#!/usr/bin/python
import spiceypy as spice
import numpy as np
import time
import sys
import os


class juno_kernel_loader:

    def __init__(self, dir, timestr='', verbose = True) : 
        self.verbose = verbose
        self.ckdir = dir + '/ck/'
        self.spkdir = dir + '/spk/'
        self.kernels = ["/fk/juno_v09.tf", "/ik/juno_jiram_v01.ti", "/ik/juno_struct_v01.ti", 
                     "/lsk/naif0012.tls", "/pck/pck00010.tpc", "/sclk/JNO_SCLKSCET.00058.tsc", 
                     "/spk/de436s.bsp", "/spk/juno_struct_v02.bsp", "/spk/jup310.bsp"]        
        for k in self.kernels : spice.furnsh(dir + k)
        if (timestr != '') : self.load(timestr)
        
    def load(self, timestr) :
        self.et = spice.str2et(str(timestr))
        self.load_ck()
        self.load_spk()        

    def et_from_YYMMDD(self, a, hms='00:00:00') : 
        timestr = str(a[0:2] + '-' + a[2:4] + '-' + a[4:6])
        if (int(a[0:2]) > 50) : timestr = '19' + timestr
        else : timestr = '20' + timestr
        et = spice.str2et(timestr + ' ' + hms)
        return et

    def load_ck(self) :     
        files = sorted(os.listdir(self.ckdir))
        self.ckfiles = []
        for file in files : 
            if (('_rec_' in file) and (not 'lbl' in file)) : 
                parts = file.split('_')
                if (len(parts) < 6) : continue
                if ((len(parts[3]) != 6) or (len(parts[4]) != 6)) : continue
                etstart = self.et_from_YYMMDD(parts[3]) #- 3600.0 * 24
                etstop  = self.et_from_YYMMDD(parts[4], hms = '23:59:59')# + 3600.0 * 24
                if ((etstart <= self.et) and (etstop >= self.et)) : 
                    if (self.verbose) : print('Loading CK file:  ' + file)
                    spice.furnsh(self.ckdir + file)
                    self.kernels.append(file)
                
    def load_spk(self) :     
        files = sorted(os.listdir(self.spkdir))
        self.ckfiles = []
        for file in files : 
            if (('_rec_' in file) and (not 'lbl' in file)) : 
                parts = file.split('_')
                if (len(parts) < 5) : continue
                if ((len(parts[2]) != 6) or (len(parts[3]) != 6)) : continue
                etstart = self.et_from_YYMMDD(parts[2])# - 3600.0 * 24
                etstop  = self.et_from_YYMMDD(parts[3], hms = '23:59:59') #+ 3600.0 * 24
                if ((etstart <= self.et) and (etstop >= self.et)) : 
                    if (self.verbose) : print('Loading SPK file: ' + file)
                    spice.furnsh(self.spkdir + file)
                    self.kernels.append(file)

class instrument_engine : 
    def __init__(self) :

        self.emission_altitude  = 0.0
        self.time_offset        = 0.0
        self.iref               = 'J2000'
        self.abcorr             = 'LT'        
        self.data               = {}
        self.constants()
        
        self.sc_id              = spice.bodn2c(self.spacecraft)
        self.target_id          = spice.bodn2c(self.target)
        self.framestring        = 'IAU_' + self.target

        self.set_emission_altitude(self.emission_altitude)

        # The angle from the centre of the pixel to the edge
        self.xpxoffset = 0.5 * self.pxscale / 1000.0
        self.ypxoffset = 0.5 * self.pyscale / 1000.0

        # Generate a key mapping, so that we know which index is what. 
        self.keys = ['lat', 'lon', 'limb_lat', 'limb_lon', 'lat_graphic', 'lon_graphic','phase', 'emission', 'incidence', 'localtime', 'losdist', 'limb_distance', 'bodyintercept']

        self.map = {}
        for key, value in enumerate(self.keys) : self.map[value] = key

        self.angles = ['lat', 'lon', 'limb_lat', 'limb_lon', 'lat_graphic', 'lon_graphic','phase', 'emission', 'incidence']


    def set_emission_altitude(self, emission_altitude) : 
        """Set the altitude of the reference spheroid, relative to IAU surface."""        

        self.emission_altitude = emission_altitude
        # Get the radius of the planet + optional altitude offset
        self.radii = spice.bodvar(self.target_id, 'RADII', 3)
        self.radii[0] = self.radii[0] + self.emission_altitude * self.radii[1] / self.radii[2]
        self.radii[1] = self.radii[1] + self.emission_altitude * self.radii[1] / self.radii[2]
        self.radii[2] = self.radii[2] + self.emission_altitude
        self.flattening = ( self.radii[0] - self.radii[2] ) / self.radii[0]

    def sc_pos(self, timestr):

        # Convert the timestring to ephemeris time
        self.et = spice.str2et(str(timestr))
        self.et = self.et + self.time_offset

        self.data['sc_wrt_planet'], self.data['target_light_time'] = spice.spkpos(self.spacecraft, self.et, self.iref, 'NONE', self.target)
        self.data['target_distance'] = spice.vnorm(self.data['sc_wrt_planet'])
        self.data['sc_wrt_sun'], light_times = spice.spkpos(self.spacecraft, self.et, self.iref, 'NONE', 'SUN')

        #self.ticks = spice.sce2t(self.sc_id, self.et)

        # Retrieve the position of the spacecraft relative to the target 
        state, self.ltcorr = spice.spkezr(self.target, self.et, self.iref, self.abcorr, self.spacecraft)
        scloc = state[0:3]

        # Get the sub-spacecraft coordinates
        point, epoch, vector = spice.subpnt('NEAR POINT/ELLIPSOID', self.target, self.et, self.framestring, self.abcorr, self.spacecraft)
        self.distance, self.data['lon_sc'], self.data['lat_sc'] = spice.reclat(point)

        # Get the localtime of the spaceship        
        hr, min, sc, time, ampm = spice.et2lst(self.et, self.target_id, self.data['lon_sc'], 'PLANETOCENTRIC')
        self.data['localtime_sc'] = hr + min / 60.0 + sc / 3600.0
        
        # Get the subsolar coordinates
        point, epoch, vector = spice.subslr('NEAR POINT/ELLIPSOID', self.target, self.et, self.framestring, self.abcorr, self.spacecraft)
        distance, self.data['lon_sun'], self.data['lat_sun'] = spice.reclat(point)

        # Convert angles from radians to degrees
        for key in ['lat_sc', 'lon_sc', 'lat_sun', 'lon_sun'] : self.data[key] = np.rad2deg(self.data[key])
        for key in ['lon_sc', 'lon_sun'] : self.data[key] = self.data[key] % 360
	
        state_planet = spice.spkssb(self.target_id, self.et, self.iref)
        inert2p = spice.pxform(self.iref, self.framestring, self.et)
        self.scloc = np.matmul(-inert2p, scloc) 

        self.i2p = spice.pxform(self.instrument, self.framestring, self.et)
        
        return self.data

    def make_pixel_vector(self, x, y, pos=0) :     

        # Define some unit vectors
        xm = np.array([1.0, 0.0, 0.0])
        ym = np.array([0.0, 1.0, 0.0])
        zm = np.array([0.0, 0.0, 1.0])

        xangle = -1.0 * self.pxscale * (x - self.middle_x + 0.5) / 1000.0
        yangle = 1.0 * self.pyscale * (y - self.middle_y + 0.5) / 1000.0
                
        rotation_x = spice.axisar(xm, xangle)
        rotation_y = spice.axisar(ym, yangle)
        i2px = np.dot(rotation_x.T, rotation_y.T).T

        # This is the central vector of the pixel
        boresight = np.dot(i2px.T, zm.T).T

        # Generate vectors for the centre and each corner
        if (pos == 0) :
            vector = boresight
        elif (pos == 1) : 
            w_pos_mat = spice.axisar(xm, self.xpxoffset)
            h_pos_mat = spice.axisar(ym, self.ypxoffset)
            corner_rot = np.matmul(w_pos_mat, h_pos_mat)
        elif (pos == 2) : 
            w_neg_mat = spice.axisar(xm, -1.0*self.xpxoffset)
            h_pos_mat = spice.axisar(ym, self.ypxoffset)
            corner_rot = np.matmul(w_neg_mat, h_pos_mat)
        elif (pos == 3) : 
            w_neg_mat = spice.axisar(xm, -1.0*self.xpxoffset)
            h_neg_mat = spice.axisar(ym, -1.0*self.ypxoffset)
            corner_rot = np.matmul(w_neg_mat, h_neg_mat)
        elif (pos == 4) : 
            w_pos_mat = spice.axisar(xm, self.xpxoffset)
            h_neg_mat = spice.axisar(ym, -1.0*self.ypxoffset)
            corner_rot = np.matmul(w_pos_mat, h_neg_mat)

        if (pos > 0) : vector = np.dot(corner_rot.T, boresight.T).T
        
        # Rotate from the centre
        vector = np.matmul(self.i2p, vector)

        return vector

    def make_pixel(self, boresights) : 
        print("None")
        
    def extract_param(self, key, pixel = 0) :
#        index = self.keys.index(key)
        index = self.map[key]
        if (pixel == 0) : arr = self.centres[:,:,index]
        if (pixel == 1) : arr = self.corners[:self.ny,:self.nx,index]
        if (pixel == 2) : arr = self.corners[1:,:self.nx,index]
        if (pixel == 3) : arr = self.corners[1:,1:,index]
        if (pixel == 4) : arr = self.corners[:self.ny,1:,index]
        return arr

    def nbr_vector_params(self) : 
        vector = self.make_pixel_vector(0, 0)
        values = self.vector_params(vector)
        return values.size
        
    def get_centres(self) :
    
        self.centres = np.zeros(shape=[self.nx, self.ny, self.nbr_vector_params()])
        self.iter = np.zeros(shape=[self.nx, self.ny])
        for (px, py), value in np.ndenumerate(self.iter):        
            vector = self.make_pixel_vector(px, py)
            self.centres[px][py] = self.vector_params(vector)

		# Calculate the max/min/avg resolution
        a = self.calculate_body_resolutions(self.centres)

        return self.centres

    def get_corners(self) : 
        """
        
        
        """
        self.corners = np.zeros(shape=[self.ny + 1, self.nx + 1, self.nbr_vector_params()])
        self.iter = np.zeros(shape=[self.ny + 1, self.nx + 1])
        for (py, px), value in np.ndenumerate(self.iter):        
            vector = self.make_pixel_vector(px, py, pos = 1)
            self.corners[py][px] = self.vector_params(vector)

		# Calculate the max/min/avg resolution
        a = self.calculate_body_resolutions(self.corners)

        return self.corners

    def calculate_body_resolutions(self, data) : 
        self.resolutions = np.zeros(shape=[data.shape[0] - 1, data.shape[1] - 1])
        for (py, px), value in np.ndenumerate(self.resolutions): 
            lon1 = np.deg2rad(data[py, px, self.map['lon']])
            lat1 = np.deg2rad(data[py, px, self.map['lat']])
            lon2 = np.deg2rad(data[py + 1, px + 1, self.map['lon']])
            lat2 = np.deg2rad(data[py + 1, px + 1, self.map['lat']])
            if (np.deg2rad(-1000.0) in [lon1, lat1, lon2, lat2]) : continue
            point1 = spice.srfrec(self.target_id, lon1, lat1)
            point2 = spice.srfrec(self.target_id, lon2, lat2)
            distance = spice.vdist(point1, point2)
            if (distance) : self.resolutions[py, px] = distance
        # This may not work
        try: 
            ontarget = self.resolutions[np.nonzero(self.resolutions)]
        except: 
            ontarget = [-1000.0]
        # if not ontarget: ontarget = [-1000.0]
        self.data['res_max'] = np.max(ontarget)
        self.data['res_min'] = np.min(ontarget)
        self.data['res_avg'] = np.mean(ontarget)
        print(self.data['res_avg'])
        return self.resolutions

    def vector_params(self, vector) : 

        ret = np.full(shape = [len(self.keys)], fill_value = -1000.0)
        
        # Calculate the distance of the boresight and the center of the target
        origin = np.array([0.0, 0.0, 0.0])
        nearpoint, rayradius = spice.nplnpt(self.scloc, vector, origin)

        # Get the intercept to calculate the radius of of the target 
        normal = spice.surfpt(origin, 1000.0 * nearpoint, self.radii[0], self.radii[1], self.radii[2])
        ret[self.map['limb_distance']] = rayradius - spice.vnorm(normal)
    
        # Calculate the 'occultation' latitude and longitude 
        ret[self.map['bodyintercept']], ret[self.map['limb_lat']], ret[self.map['limb_lon']] = spice.reclat(nearpoint) 

        # Test if the boresight intersects with our target surface 
        try: 
            point = spice.surfpt(self.scloc, vector, self.radii[0], self.radii[1], self.radii[2])
            intercept = True
        except: 
            intercept = False
    
        if (intercept) : 
        
            # Get some angles 
            ret[self.map['phase']], ret[self.map['incidence']], ret[self.map['emission']] = spice.illum(self.target, self.et , self.abcorr, self.spacecraft, point)

            # Calculate the planetocentric coordinates 
            distance, ret[self.map['lon']], ret[self.map['lat']] = spice.reclat(point)
            ret[self.map['losdist']] = spice.vnorm(self.scloc - point)
            
             # Calculate the planetographic coordinates
            ret[self.map['lon_graphic']], ret[self.map['lat_graphic']], bodyintercept = spice.recpgr(self.target, point, self.radii[0], self.flattening)
            
            # Get the localtime, and convert to decimal hours 
            hr, min, sc, time, ampm = spice.et2lst(self.et, self.target_id, ret[self.map['lon']], 'PLANETOCENTRIC')
            ret[self.map['localtime']] = hr + min / 60.0 + sc / 3600.0
            
        # For the angles, convert radians to degrees
        for key in self.angles : 
            if (ret[self.map[key]] != -1000.0) : 
                    ret[self.map[key]] = np.rad2deg(ret[self.map[key]])

        # Makes sure longitudes wrap 0 to 360, spice returns the Earth-like -180 to 180. 
        longitudes = ['lon', 'limb_lon', 'lon_graphic']
        for key in longitudes : 
            ret[self.map[key]] = ret[self.map[key]] % 360

        return np.array(ret)
                   
    def constants(self) : 
        # This function should be overridden to define the geometric parameters of the instrument.
        self.spacecraft = 'JUNO'
        
        
class juno_jiram_geometry(instrument_engine) : 
    def constants(self) :     

        # Who are we and what are we looking at
        self.spacecraft     = 'JUNO'
        self.target         = 'JUPITER'
        self.instrument     = 'JUNO_JIRAM_I_MBAND'

        # Number of pixels in the x and y directions
        self.nx             = 432
        self.ny             = 128

        # The angular size of each pixel (mrad)
        self.pxscale        = 0.24
        self.pyscale        = 0.24 
        
        # The pixel that defines  the centre of the instrument boresight
        self.middle_y       = 128.0 / 2.0 #(128.0) / 2.0
        self.middle_x       = (432.0) / 2.0
        
    def set_mode(self, mode_id, imager = True) :
        """ 
            JIRAM has three imaging modes (L, M, and L + M) and three 
            spectral modes (16, 64, and 256 spatial pixels). Only the 
            full spectral mode with 256 pixels is implemented here.
            
            `mode_id` is the `INSTRUMENT_MODE_ID` string in the LBL metadata file. 
            It can look something like this: SCI_I2_S3
        """
        if ('I1' in mode_id) : 
            self.instrument = 'JUNO_JIRAM_I'
            self.nx = 432
            self.ny = 266
            self.middle_x = self.nx / 2.0
            self.middle_y = self.ny / 2.0
        elif ('I2' in mode_id) : self.instrument = 'JUNO_JIRAM_I_MBAND'
        elif ('I3' in mode_id) : self.instrument = 'JUNO_JIRAM_I_LBAND'
        if (imager == False) : 
            self.instrument = 'JUNO_JIRAM_S'
            self.nx = 256
            self.ny = 1
            self.middle_y = 0.5
            self.middle_x = 256.0 / 2.0 + 1
            

#=========================================================================================
#=========================================================================================
#=========================================================================================
#=========================================================================================
#
# Below is a slightly modified version of https://github.com/mkelley/pds3, required to
# read the JIRAM images stored on the PDS. 
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = [
    'read_label',
    'read_ascii_table',
    'read_image',
    'read_table'
]

from collections import OrderedDict

try:
    from ply import lex, yacc
except ImportError:
    raise ImportError("pds3 requires the PLY (Python Lex-Yacc) module.")

class IllegalCharacter(Exception):
    pass

class PDS3Keyword(str):
    """PDS3 keyword.

    In the following, the keyword is "IMAGE":
      
      OBJECT = IMAGE
        ...
      END_OBJECT = IMAGE

    """
    def __new__(cls, value):
        return str.__new__(cls, value)

class PDS3Object(OrderedDict):
    """PDS3 data object definition.

    OBJECT = IMAGE
      ...
    END_OBJECT = IMAGE

    """
    pass
    
class PDS3Group(OrderedDict):
    """PDS3 group statement.

    GROUP = SHUTTER_TIMES
      ...
    END_GROUP = SHUTTER_TIMES

    """
    pass
    

PDS3_DATA_TYPE_TO_DTYPE = {
    'IEEE_REAL': '>f',
    'LSB_INTEGER': '<i',
    'LSB_UNSIGNED_INTEGER': '<u',
    'MAC_INTEGER': '>i',
    'MAC_REAL': '>f',
    'MAC_UNSIGNED_INTEGER': '>u',
    'MSB_UNSIGNED_INTEGER': '>u',
    'MSB_INTEGER': '>i',
    'PC_INTEGER': '<i',
    'PC_UNSIGNED_INTEGER': '<u',
    'SUN_INTEGER': '>i',
    'SUN_REAL': '>f',
    'SUN_UNSIGNED_INTEGER': '>u',
    'VAX_INTEGER': '<i',
    'VAX_UNSIGNED_INTEGER': '<u',
}

class PDS3Parser():
    tokens = ['KEYWORD', 'POINTER', 'STRING', 'INT', 'REAL',
              'UNIT', 'DATE', 'END']

    literals = list('=(){},')

    t_POINTER = r'\^[A-Z0-9_]+'
    t_ignore_COMMENT = r'/\*.+?\*/'
    t_ignore = ' \t\r\n'

    # lower case PDS3 to astropy unit translation
    unit_translate = dict(v='V', k='K')

    def __init__(self, debug=False):
        import os
        self.debug = debug
        self.lexer = lex.lex(module=self, debug=self.debug)
        self.parser = yacc.yacc(module=self, debug=self.debug,
                                write_tables=0)

    def parse(self, raw_label, **kwargs):
        return self.parser.parse(raw_label, lexer=self.lexer, debug=self.debug,
                                 **kwargs)

    def t_KEYWORD(self, t):
        r'[A-Z][A-Z0-9_:]+'
        if t.value == 'END':
            t.type = 'END'
        return t

    def t_DATE(self, t):
        r'\d\d\d\d-\d\d-\d\d(T\d\d:\d\d(:\d\d(.\d+)?)?)?Z?'
        from astropy.time import Time
        t.value = Time(t.value, scale='utc')
        return t

    def t_UNIT(self, t):
        r'<[\w*^\-/]+>'
        import astropy.units as u

        # most astropy units are lower-case versions of the PDS3 units
        unit = t.value[1:-1].lower()

        # but not all
        if unit in self.unit_translate:
            unit = self.unit_translate[unit]

        t.value = u.Unit(unit)
        return t

    def t_STRING(self, t):
        r'"[^"]+"'
        t.value = t.value[1:-1].replace('\r', '')
        return t
 
    def t_REAL(self, t):
        r'[+-]?(([0-9]+\.[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?'
        t.value = float(t.value)
        return t

    def t_INT(self, t):
        r'[+-]?[0-9]+'
        t.value = int(t.value)
        return t

    def t_error(self, t):
        raise IllegalCharacter(t.value[0])

    def lexer_test(self, data):
        self.lexer.input(data)
        while True:
            tok = self.lexer.token()
            if not tok:
                break
            print(tok)

    def p_label(self, p):
        """label : record
                 | label record
                 | label END"""
        if len(p) == 2:
            # record
            p[0] = [p[1]]
        elif p[2] == 'END':
            # label END
            p[0] = p[1]
        else:
            # label record
            p[0] = p[1] + [p[2]]

    def p_record(self, p):
        """record : KEYWORD '=' value
                  | POINTER '=' INT
                  | POINTER '=' STRING
                  | POINTER '=' '(' STRING ',' INT ')'"""
        if len(p) == 4:
            p[0] = (p[1], p[3])
        else:
            p[0] = (p[1], (p[4], p[6]))

    def p_value(self, p):
        """value : STRING
                 | DATE
                 | KEYWORD
                 | number
                 | pds_set
                 | quantity
                 | sequence"""
        p[0] = p[1]

    def p_value_quantity(self, p):
        """quantity : number UNIT"""
        p[0] = p[1] * p[2]

    def p_number(self, p):
        """number : INT
                  | REAL"""
        p[0] = p[1]

    def p_pds_set(self, p):
        """pds_set : '{' value '}'
                   | '{' sequence_values '}'"""
        p[0] = set(p[2])

    def p_sequence(self, p):
        """sequence : '(' value ')'
                    | '(' sequence_values ')'
                    | '{' value '}'
                    | '{' sequence_values '}'"""
        p[0] = p[2]

    def p_sequence_values(self, p):
        """sequence_values : value ','
                           | sequence_values value ','
                           | sequence_values value"""
        if p[2] == ',':
            p[0] = [p[1]]
        else:
            p[0] = p[1] + [p[2]]

    def p_error(self, p):
        if p:
            print("Syntax error at '%s'" % p.value)
        else:
            print("Syntax error at EOF")


def _records2dict(records, object_index=0):
    """Convert a list of PDS3 records to a dictionary.

    Parameters
    ----------
    records : list
      List of key-item pairs.
    object_index : int, optional
      Extract just a single object or group, starting at the index
      `object_index`.

    Returns
    -------
    d : dict
      The dictionary.
    last_index : int, optional
      When extracting a single object or group, also return the last
      index of that object.

    """
    
    from collections import OrderedDict
    
    label = OrderedDict()
    start = 0
    if object_index != 0:
        start = object_index
        object_name = records[start][1]
        start += 1
    else:
        object_name = None

    i = start
    while i < len(records):
        # groups and objects are both terminated with 'END_...'
        if (records[i] == ('END_OBJECT', object_name)
            or records[i] == ('END_GROUP', object_name)):
            return label, i
        elif records[i][0] in ['OBJECT', 'GROUP']:
            key = PDS3Keyword(records[i][1])
            value, j = _records2dict(records, i)
            if records[i][0] == 'OBJECT':
                value = PDS3Object(value)
            elif records[i][0] == 'GROUP':
                value = PDS3Group(value)
            i = j
        else:
            key = records[i][0]
            value = records[i][1]

        if key in label:
            if not isinstance(label[key], list):
                label[key] = [label[key]]
            label[key].append(value)
        else:
            label[key] = value
        i += 1

    return label

def _find_file(filename, path='.'):
    return filename
    """Search a directory for a file.

    PDS3 file names are required to be upper case, but in case-
    sensitive file systems, the PDS3 file name many not match the file
    system file name.  `_find_file` will make a case insentive
    comparison to all files in the given path.

    Parameters
    ----------
    filename : string
      The PDS3 file name for which to search.
    path : string, optional
      Search within this directory path.

    Returns
    -------
    fn : string
      The actual name of the file on the file system, including the
      path.

    Raises
    ------
    IOError
      When the file is not found.
      
    """

    import os
    for f in sorted(os.listdir(path)):
        if f.lower() == filename.lower():
            f = os.path.sep.join([path, f])
            return f
    raise IOError("No match for {} in {}".format(filename, path))

def read_label(filename, debug=False):
    """Read in a PDS3 label.

    Parameters
    ----------
    filename : string
      The name of the file to read.

    Returns
    -------
    label : dict
      The label as a `dict`.

    Raises
    ------
    IllegalCharacter

    Notes
    -----
    Objects and groups are returned as dictionaries containing all
    their sub-keywords.  Multiple objects (e.g., columns) with the
    same name are returned as a `list` of objects.

    """

    raw_label = ''
    with open(filename, 'rb') as inf:
        while True:
            line = inf.readline()
            raw_label += line.decode('ascii')
            if line.strip() == b'END' or line == '':
                break

    parser = PDS3Parser(debug=debug)
    records = parser.parse(raw_label)
    return _records2dict(records)

def read_ascii_table(label, key, path='.'):
    """Read an ASCII table as described by the label.

    Only fixed length records are supported.

    Parameters
    ----------
    label : dict
      The label, as read by `read_label`.
    key : string
      The label key of the object that describes the table.
    path : string, optional
      Directory path to label/table.

    Returns
    -------
    table : astropy Table

    Raises
    ------
    NotImpementedError
    ValueError

    """

    import numpy as np
    from astropy.io import ascii

    # The table object description.
    desc = label[key]

    if not isinstance(desc['COLUMN'], list):
        # For tables with a single column, desc['COLUMN'] needs to be a list
        desc['COLUMN'] = [desc['COLUMN']]

    # Setup table column formats
    n = desc['COLUMNS']
    col_starts = []
    col_ends = []
    converters = dict()
    def repeat(dtype, col):
        n = col.get('ITEMS', 1)
        return dtype if n == 1 else (dtype,) * n
    
    for i in range(n):
        col = desc['COLUMN'][i]
        col_starts.append(col['START_BYTE'] - 1)
        col_ends.append(col_starts[-1] + col['BYTES'] - 1)

        if col['DATA_TYPE'] == 'ASCII_REAL':
            dtype = repeat(np.float, col)
        elif col['DATA_TYPE'] == 'ASCII_INTEGER':
            dtype = repeat(np.int, col)
        elif col['DATA_TYPE'] == 'CHARACTER':
            dtype = repeat('S{}'.format(col['BYTES']), col)
        else:
            raise ValueError("Unknown data type: ", col['DATA_TYPE'])
        converters['col{}'.format(i+1)] = [ascii.convert_numpy(dtype)]

    nrows = desc['ROWS']

    # Open the file object, and skip ahead to the start of the table,
    # if needed.  Read the table.
    if isinstance(label['^' + key], tuple):
        filename, start = label['^' + key]
        start = int(start) - 1
        filename = _find_file(filename, path=path)
        if 'RECORD_BYTES' in label:
            record_bytes = label['RECORD_BYTES']
        else:
            record_bytes = desc['RECORD_BYTES']

        #inf = open(filename, 'r')
        #inf.seek(record_bytes * start)
    else:
        filename = _find_file(label['^' + key], path=path)
        start = 0
        #inf = open(filename, 'r')

    table = ascii.read(filename, format='fixed_width_no_header',
                       data_start=start, data_end=nrows+start,
                       col_starts=col_starts, col_ends=col_ends,
                       converters=converters, guess=False)
    #inf.close()

    # Mask data
    for i in range(n):
        col = desc['COLUMN'][i]
        missing_constant = col.get('MISSING_CONSTANT', None)
        if missing_constant is None:
            continue

        j = table.columns[i] == missing_constant
        if np.any(j):
            table.columns[i].mask = j

    # Save column meta data.
    for i in range(n):
        for j in range(desc.get('ITEMS', 1)):
            col = desc['COLUMN'][i]
            table.columns[i].name = col['NAME']
            if 'DESCRIPTION' in col:
                table.columns[i].description = col['DESCRIPTION']

    # Save table meta data.
    for k, v in desc.items():
        if k is not 'COLUMN':
            table.meta[k] = v

    return table

def read_binary_table(label, key, path='.'):
    """Read a binary table as described by the label.

    Parameters
    ----------
    label : dict
      The label, as read by `read_label`.
    key : string
      The label key of the object that describes the table.
    path : string, optional
      Directory path to label/table.

    Returns
    -------
    table : 

    Raises
    ------
    NotImpementedError
    ValueError

    """

    import numpy as np

    desc = label[key]
    dtype = []
    byte = 1
    offset = np.zeros(len(desc['COLUMN']))
    scale = np.ones(len(desc['COLUMN']))
    for i, col in enumerate(desc['COLUMN']):
        if col['START_BYTE'] != byte:
            raise NotImplementedError("Table requires skipping bytes.")

        if col['ITEMS'] != 1:
            raise NotImplementedError("Table requires mutliple items per column.")

        byte += col['START_BYTE']
        x = PDS3_DATA_TYPE_TO_DTYPE[col['DATA_TYPE']] + str(col['ITEM_BYTES'])
        dtype.append((col['NAME'], x))

        offset[i] = col.get('OFFSET', 0.0)
        scale[i] = col.get('SCALE', 1.0)

    if desc['ROW_BYTES'] != bytes:
        raise NotImplementedError("Table requires skipping bytes.")

    dtype = np.dtype(dtype)

    if isinstance(label['^' + key], tuple):
        filename, start = label['^' + key]
        start = int(start) - 1
        filename = _find_file(filename, path=path)
        if 'RECORD_BYTES' in label:
            record_bytes = label['RECORD_BYTES']
        else:
            record_bytes = desc['RECORD_BYTES']

        inf = open(filename, 'r')
        inf.seek(record_bytes * start)
    else:
        filename = _find_file(label['^' + key], path=path)
        start = 0
        inf = open(filename, 'r')

    n = desc['ROWS']
    data = np.fromfile(inf, dtype=dtype, count=n).view(np.recarray)
    return data

def read_image(filename, label, key, path=".", scale_and_offset=True, verbose=False):
    """Read an image as described by the label.

    The image is not reordered for display orientation.

    When there are interleaved data (i.e., the BANDS keyword is
    present), they will be separated and a data cube returned.

    Parameters
    ----------
    label : dict
      The label, as read by `read_label`.
    key : string
      The label key of the object that describes the image.
    path : string, optional
      Directory path to label/table.
    scale_and_offset : bool, optional
      Set to `True` to apply the scale and offset factors and return
      floating point data.
    verbose : bool, optional
      Print some informational info.

    Returns
    -------
    im : ndarray
      The image.  If there are multiple bands, the first index
      iterates over each band.

    """

    import os.path
    import warnings
    import numpy as np

#    warnings.warn("read_image is a basic and incomplete reader.")

    # The image object description.
    desc = label['FILE']['IMAGE']

    #print(desc)

    shape = np.array((desc['LINES'], desc['LINE_SAMPLES']))
    size = desc['SAMPLE_BITS'] // 8

    if 'BANDS' in desc:
        bands = desc['BANDS']
    else:
        bands = 1

    if 'LINE_PREFIX_BYTES' in desc:
        prefix_shape = (shape[0], desc['LINE_PREFIX_BYTES'])
    else:
        prefix_shape = (0, 0)

    if 'LINE_SUFFIX_BYTES' in desc:
        suffix_shape = (shape[0], desc['LINE_SUFFIX_BYTES'])
    else:
        suffix_shape = (0, 0)

    line_size = prefix_shape[0] + shape[0] * size + suffix_shape[0]
    #print(PDS3_DATA_TYPE_TO_DTYPE)
    if desc['SAMPLE_TYPE'] in PDS3_DATA_TYPE_TO_DTYPE:
        dtype = '{}{:d}'.format(PDS3_DATA_TYPE_TO_DTYPE[desc['SAMPLE_TYPE']],
                                size)
    else:
        raise NotImplemented('SAMPLE_TYPE={}'.format(desc['SAMPLE_TYPE']))
   # filename = label['FILE']['^IMAGE']
    dtype = 'f4'
    start_record = 1
#    filename, start_record = label['^{}'.format(key)]
    #print(filename, start_record)
    found_filename = _find_file(filename)
    start = (start_record - 1) * int(label['FILE']['RECORD_BYTES'])
    verbose = 1
   # if verbose:
        #print('''Image shape: {} samples
#Number of bands: {}
#Stored line size, including prefix, and suffix: {} bytes
#Numpy data type: {}
#Filename: {} ({})
#Start byte: {}'''.format(shape, line_size, bands, dtype, filename,
#                         found_filename, start))
    # The file is read into a an array of bytes in order to properly
    # handle line prefix and suffix.
    with open(found_filename, 'rb') as inf:
        inf.seek(0)
        n = np.prod(line_size * shape[1])
        buf = inf.read(n)
        if n != len(buf):
            raise IOError("Expected {} bytes of data, but only read {}".format(
                    n, len(buf)))
    
    # remove prefix and suffix, convert to image data type and shape
    data = np.frombuffer(buf, dtype=np.uint8, count=n)
    del buf
    data = data.reshape((line_size, shape[1]))
    s = slice(prefix_shape[0], (suffix_shape[0] if suffix_shape[0] > 0 else None))
    data = data[s, :].flatten()
    im = np.frombuffer(data, dtype=dtype, count=np.prod(shape)).reshape(shape)
    del data

    # separate out interleaved data
    if bands > 1:
        if desc['BAND_STORAGE_TYPE'].upper() == 'SAMPLE_INTERLEAVED':
            im = im.reshape((shape[0], shape[1] // bands, bands))
            im = np.rollaxis(im, 2)
        elif desc['BAND_STORAGE_TYPE'].upper() == 'LINE_INTERLEAVED':
            im = im.reshape((shape[0] // bands, bands, shape[1]))
            im = np.rollaxis(im, 1)
        else:
            raise ValueError('Incorrect BAND_STORAGE_TYPE: {}'.format(
                    desc['BAND_STORAGE_TYPE']))

    if ('OFFSET' in desc or 'SCALING_FACTOR' in desc) and scale_and_offset:
        im = (desc.get('OFFSET', 0)
              + im.astype(float) * desc.get('SCALING_FACTOR', 1))

    return im

def read_table(label, key, path='.'):
    """Read table as described by the label.

    Calls `read_ascii_table` or `read_binary_table` as appropriate.

    Parameters
    ----------
    label : dict
      The label, as read by `read_label`.
    key : string
      The label key of the object that describes the table.
    path : string, optional
      Directory path to label/table.

    Returns
    -------
    table : astropy Table

    """

    format = label[key]['INTERCHANGE_FORMAT']
    if format == 'ASCII':
        return read_ascii_table(label, key, path=path)
    else:
        raise NotImplementedError("Table format not implemented: {}".format(format))
