import re
from datetime import datetime
from math import pi, pow
from sgp4.ext import days2mdhms, jday
from sgp4.model import Satellite
from sgp4.propagation import sgp4init

INT_RE = re.compile(r'[+-]?\d*')
FLOAT_RE = re.compile(r'[+-]?\d*(\.\d*)?')

LINE1 = '1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN'
LINE2 = '2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN'

def twoline2rv(longstr1, longstr2, whichconst, afspc_mode=False):
    deg2rad  =   pi / 180.0;         #    0.0174532925199433
    xpdotp   =  1440.0 / (2.0 *pi);  #  229.1831180523293

    tumin = whichconst.tumin

    satrec = Satellite()
    satrec.error = 0;
    satrec.whichconst = whichconst  # Python extension: remembers its consts

    line = longstr1.rstrip()
    # try/except is not well supported by Numba
    if (len(line) >= 64 and
        line.startswith('1 ') and
        line[8] == ' ' and
        line[23] == '.' and
        line[32] == ' ' and
        line[34] == '.' and
        line[43] == ' ' and
        line[52] == ' ' and
        line[61] == ' ' and
        line[63] == ' '):

        _saved_satnum = satrec.satnum = int(line[2:7])
        # classification = line[7] or 'U'
        # intldesg = line[9:17]
        two_digit_year = int(line[18:20])
        satrec.epochdays = float(line[20:32])
        satrec.ndot = float(line[33:43])
        satrec.nddot = float(line[44] + '.' + line[45:50])
        nexp = int(line[50:52])
        satrec.bstar = float(line[53] + '.' + line[54:59])
        ibexp = int(line[59:61])
        # numb = int(line[62])
        # elnum = int(line[64:68])
    else:
        raise ValueError(error_message.format(1, LINE1, line))

    line = longstr2.rstrip()
    if (len(line) >= 69 and
        line.startswith('2 ') and
        line[7] == ' ' and
        line[11] == '.' and
        line[16] == ' ' and
        line[20] == '.' and
        line[25] == ' ' and
        line[33] == ' ' and
        line[37] == '.' and
        line[42] == ' ' and
        line[46] == '.' and
        line[51] == ' '):

        satrec.satnum = int(line[2:7])
        if _saved_satnum != satrec.satnum:
            raise ValueError('Object numbers in lines 1 and 2 do not match')

        satrec.inclo = float(line[8:16])
        satrec.nodeo = float(line[17:25])
        satrec.ecco = float('0.' + line[26:33].replace(' ', '0'))
        satrec.argpo = float(line[34:42])
        satrec.mo = float(line[43:51])
        satrec.no = float(line[52:63])
        #revnum = line[63:68]
    #except (AssertionError, IndexError, ValueError):
    else:
        raise ValueError(error_message.format(2, LINE2, line))

    #  ---- find no, ndot, nddot ----
    satrec.no   = satrec.no / xpdotp; #   rad/min
    satrec.nddot= satrec.nddot * pow(10.0, nexp);
    satrec.bstar= satrec.bstar * pow(10.0, ibexp);

    #  ---- convert to sgp4 units ----
    satrec.a    = pow( satrec.no*tumin , (-2.0/3.0) );
    satrec.ndot = satrec.ndot  / (xpdotp*1440.0);  #   ? * minperday
    satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);

    #  ---- find standard orbital elements ----
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo  * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     * deg2rad;

    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;
	
    if two_digit_year < 57:
        year = two_digit_year + 2000;
    else:
        year = two_digit_year + 1900;

    mon,day,hr,minute,sec = days2mdhms(year, satrec.epochdays);
    sec_whole, sec_fraction = divmod(sec, 1.0)

    satrec.epochyr = year
    satrec.jdsatepoch = jday(year,mon,day,hr,minute,sec);
    satrec.epoch = datetime(year, mon, day, hr, minute, int(sec_whole),
                            int(sec_fraction * 1000000.0 // 1.0))

    #  ---------------- initialize the orbit at sgp4epoch -------------------
    sgp4init(whichconst, afspc_mode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar,
             satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no,
             satrec.nodeo, satrec)

    return satrec

def verify_checksum(*lines):
    for line in lines:
        checksum = line[68:69]
        if not checksum.isdigit():
            continue
        checksum = int(checksum)
        computed = compute_checksum(line)
        if checksum != computed:
            complaint = ('TLE line gives its checksum as {}'
                         ' but in fact tallies to {}:\n{}')
            raise ValueError(complaint.format(checksum, computed, line))

def fix_checksum(line):
    return line[:68].ljust(68) + str(compute_checksum(line))

def compute_checksum(line):
    return sum((int(c) if c.isdigit() else c == '-') for c in line[0:68]) % 10