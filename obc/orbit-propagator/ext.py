from math import (acos, asinh, atan2, copysign, cos, fabs, fmod, pi, sin, sinh, sqrt, tan)

undefined = None

def mag(x):
     return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

def cross(vec1, vec2, outvec):
     outvec[0]= vec1[1]*vec2[2] - vec1[2]*vec2[1];
     outvec[1]= vec1[2]*vec2[0] - vec1[0]*vec2[2];
     outvec[2]= vec1[0]*vec2[1] - vec1[1]*vec2[0];

def dot(x, y):
     return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);

def angle(vec1, vec2):
     small     = 0.00000001;
     undefined = 999999.1;

     magv1 = mag(vec1);
     magv2 = mag(vec2);

     if magv1*magv2 > small*small:

         temp= dot(vec1,vec2) / (magv1*magv2);
         if fabs(temp) > 1.0:
             temp = copysign(1.0, temp)
         return acos( temp );

     else:
         return undefined;

def newtonnu(ecc, nu):
     #  ---------------------  implementation   ---------------------
     e0= 999999.9;
     m = 999999.9;
     small = 0.00000001;

     #  --------------------------- circular ------------------------
     if fabs(ecc) < small:

         m = nu;
         e0= nu;

     else:
         #  ---------------------- elliptical -----------------------
         if ecc < 1.0-small:

             sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
             cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
             e0  = atan2( sine,cose );
             m   = e0 - ecc*sin(e0);

         else:
             #  -------------------- hyperbolic  --------------------
             if ecc > 1.0 + small:

                 if ecc > 1.0 and fabs(nu)+0.00001 < pi-acos(1.0 /ecc):

                     sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                     e0  = asinh( sine );
                     m   = ecc*sinh(e0) - e0;

             else:
                 #  ----------------- parabolic ---------------------
                 if fabs(nu) < 168.0*pi/180.0:

                     e0= tan( nu*0.5  );
                     m = e0 + (e0*e0*e0)/3.0;

     if ecc < 1.0:

         m = fmod( m,2.0 *pi );
         if m < 0.0:
             m = m + 2.0 *pi;
         e0 = fmod( e0,2.0 *pi );

     return e0, m

def rv2coe(r, v, mu):

     hbar = [None, None, None]
     nbar = [None, None, None]
     ebar = [None, None, None]
     typeorbit = [None, None, None];

     twopi  = 2.0 * pi;
     halfpi = 0.5 * pi;
     small  = 0.00000001;
     undefined = 999999.1;
     infinite  = 999999.9;

     #  -------------------------  implementation   -----------------
     magr = mag( r );
     magv = mag( v );

     #  ------------------  find h n and e vectors   ----------------
     cross( r,v, hbar );
     magh = mag( hbar );
     if magh > small:

         nbar[0]= -hbar[1];
         nbar[1]=  hbar[0];
         nbar[2]=   0.0;
         magn = mag( nbar );
         c1 = magv*magv - mu /magr;
         rdotv = dot( r,v );
         for i in range(0, 3):
             ebar[i]= (c1*r[i] - rdotv*v[i])/mu;
         ecc = mag( ebar );

         #  ------------  find a e and semi-latus rectum   ----------
         sme= ( magv*magv*0.5  ) - ( mu /magr );
         if fabs( sme ) > small:
             a= -mu  / (2.0 *sme);
         else:
             a= infinite;
         p = magh*magh/mu;

         #  -----------------  find inclination   -------------------
         hk= hbar[2]/magh;
         incl= acos( hk );

         #  --------  determine type of orbit for later use  --------
         #  ------ elliptical, parabolic, hyperbolic inclined -------
         typeorbit = 'ei'
         if ecc < small:

             #  ----------------  circular equatorial ---------------
             if  incl < small or fabs(incl-pi) < small:
                 typeorbit = 'ce'
             else:
                 #  --------------  circular inclined ---------------
                 typeorbit = 'ci'

         else:

             #  - elliptical, parabolic, hyperbolic equatorial --
             if incl < small or fabs(incl-pi) < small:
                 typeorbit = 'ee'

         #  ----------  find longitude of ascending node ------------
         if magn > small:

             temp= nbar[0] / magn;
             if fabs(temp) > 1.0:
                 temp = copysign(1.0, temp)
             omega= acos( temp );
             if nbar[1] < 0.0:
                 omega= twopi - omega;

         else:
             omega= undefined;

         #  ---------------- find argument of perigee ---------------
         if typeorbit == 'ei':

             argp = angle( nbar,ebar);
             if ebar[2] < 0.0:
                 argp= twopi - argp;

         else:
             argp= undefined;

         #  ------------  find true anomaly at epoch    -------------
         if typeorbit[0] == 'e':

             nu =  angle( ebar,r);
             if rdotv < 0.0:
                 nu= twopi - nu;

         else:
             nu= undefined;

         #  ----  find argument of latitude - circular inclined -----
         if typeorbit == 'ci':

             arglat = angle( nbar,r );
             if r[2] < 0.0:
                 arglat= twopi - arglat;
             m = arglat;

         else:
             arglat= undefined;

         #  -- find longitude of perigee - elliptical equatorial ----
         if ecc > small and typeorbit == 'ee':

             temp= ebar[0]/ecc;
             if fabs(temp) > 1.0:
                 temp = copysign(1.0, temp)
             lonper= acos( temp );
             if ebar[1] < 0.0:
                 lonper= twopi - lonper;
             if incl > halfpi:
                 lonper= twopi - lonper;

         else:
             lonper= undefined;

         #  -------- find true longitude - circular equatorial ------
         if magr > small and typeorbit == 'ce':

             temp= r[0]/magr;
             if fabs(temp) > 1.0:
                 temp = copysign(1.0, temp)
             truelon= acos( temp );
             if r[1] < 0.0:
                 truelon= twopi - truelon;
             if incl > halfpi:
                 truelon= twopi - truelon;
             m = truelon;

         else:
             truelon= undefined;

         #  ------------ find mean anomaly for all orbits -----------
         if typeorbit[0] == 'e':
             e, m = newtonnu(ecc, nu);

     else:
        p    = undefined;
        a    = undefined;
        ecc  = undefined;
        incl = undefined;
        omega= undefined;
        argp = undefined;
        nu   = undefined;
        m    = undefined;
        arglat = undefined;
        truelon= undefined;
        lonper = undefined;

     return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper

def jday(year, mon, day, hr, minute, sec):
     return (367.0 * year -
          7.0 * (year + ((mon + 9.0) // 12.0)) * 0.25 // 1.0 +
          275.0 * mon // 9.0 +
          day + 1721013.5 +
          ((sec / 60.0 + minute) / 60.0 + hr) / 24.0  #  ut in days
          #  - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
          )

def days2mdhms(year, days):
     lmonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

     dayofyr = int(days // 1.0);
     #  ----------------- find month and day of month ----------------
     if (year % 4) == 0:
       lmonth = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

     i = 1;
     inttemp = 0;
     while dayofyr > inttemp + lmonth[i-1] and i < 12:

       inttemp = inttemp + lmonth[i-1];
       i += 1;

     mon = i;
     day = dayofyr - inttemp;

     #  ----------------- find hours minutes and seconds -------------
     temp = (days - dayofyr) * 24.0;
     hr   = int(temp // 1.0);
     temp = (temp - hr) * 60.0;
     minute  = int(temp // 1.0);
     sec  = (temp - minute) * 60.0;

     return mon, day, hr, minute, sec

def invjday(jd):
     #  --------------- find year and days of the year ---------------
     temp    = jd - 2415019.5;
     tu      = temp / 365.25;
     year    = 1900 + int(tu // 1.0);
     leapyrs = int(((year - 1901) * 0.25) // 1.0);

     #  optional nudge by 8.64x10-7 sec to get even outputs
     days    = temp - ((year - 1900) * 365.0 + leapyrs) + 0.00000000001;

     #  ------------ check for case of beginning of a year -----------
     if (days < 1.0):
         year    = year - 1;
         leapyrs = int(((year - 1901) * 0.25) // 1.0);
         days    = temp - ((year - 1900) * 365.0 + leapyrs);

     #  ----------------- find remaing data  -------------------------
     mon, day, hr, minute, sec = days2mdhms(year, days);
     sec = sec - 0.00000086400;
     return year, mon, day, hr, minute, sec