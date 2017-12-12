def jit(jit_this=None, **jit_options):
  if jit_this is not None:
	   def fake_jit(*args, **kwargs):
			return jit_this(*args, **kwargs)
	   return fake_jit
  else:
	   def partial_fake_jit(jit_this, **jit_options):
			return jit(jit_this, **jit_options)
	   return partial_fake_jit

from math import atan2, cos, fabs, pi, sin, sqrt

deg2rad = pi / 180.0;
_nan = float('NaN')
false = (_nan, _nan, _nan)
true = True
twopi = 2.0 * pi

def _initl(
       satn,      whichconst,
       ecco,   epoch,  inclo,   no,
       method,
       afspc_mode,
       ):

     # sgp4fix use old way of finding gst

     #  ----------------------- earth constants ----------------------
     #  sgp4fix identify constants and allow alternate values
     tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
     x2o3   = 2.0 / 3.0;

     #  ------------- calculate auxillary epoch quantities ----------
     eccsq  = ecco * ecco;
     omeosq = 1.0 - eccsq;
     rteosq = sqrt(omeosq);
     cosio  = cos(inclo);
     cosio2 = cosio * cosio;

     #  ------------------ un-kozai the mean motion -----------------
     ak    = pow(xke / no, x2o3);
     d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
     del_  = d1 / (ak * ak);
     adel  = ak * (1.0 - del_ * del_ - del_ *
             (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
     del_  = d1/(adel * adel);
     no    = no / (1.0 + del_);

     ao    = pow(xke / no, x2o3);
     sinio = sin(inclo);
     po    = ao * omeosq;
     con42 = 1.0 - 5.0 * cosio2;
     con41 = -con42-cosio2-cosio2;
     ainv  = 1.0 / ao;
     posq  = po * po;
     rp    = ao * (1.0 - ecco);
     method = 'n';

     #  sgp4fix modern approach to finding sidereal time
     if afspc_mode:
         #  sgp4fix use old way of finding gst
         #  count integer number of days from 0 jan 1970
         ts70  = epoch - 7305.0;
         ds70 = (ts70 + 1.0e-8) // 1.0;
         tfrac = ts70 - ds70;
         #  find greenwich location at epoch
         c1    = 1.72027916940703639e-2;
         thgr70= 1.7321343856509374;
         fk5r  = 5.07551419432269442e-15;
         c1p2p = c1 + twopi;
         gsto  = (thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r) % twopi
         if gsto < 0.0:
             gsto = gsto + twopi;

     else:
        gsto = _gstime(epoch + 2433281.5);

     return (
       no,
       method,
       ainv,  ao,    con41,  con42, cosio,
       cosio2,eccsq, omeosq, posq,
       rp,    rteosq,sinio , gsto,
       )

def sgp4init(
       whichconst, afspc_mode,   satn,     epoch,
       xbstar,  xecco, xargpo,
       xinclo,  xmo,   xno,
       xnodeo,  satrec,
       ):
     
     temp4    =   1.5e-12;

     #  ----------- set all near earth variables to zero ------------
     satrec.isimp   = 0;   satrec.method = 'n'; satrec.aycof    = 0.0;
     satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
     satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
     satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
     satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
     satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
     satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
     satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot  = 0.0;
     satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf   = 0.0;
	 
     satrec.bstar   = xbstar;
     satrec.ecco    = xecco;
     satrec.argpo   = xargpo;
     satrec.inclo   = xinclo;
     satrec.mo	    = xmo;
     satrec.no	    = xno;
     satrec.nodeo   = xnodeo;

     #  sgp4fix add opsmode
     satrec.afspc_mode = afspc_mode;

     #  ------------------------ earth constants -----------------------
     #  sgp4fix identify constants and allow alternate values
     tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
     ss     = 78.0 / radiusearthkm + 1.0;
     #  sgp4fix use multiply for speed instead of pow
     qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
     qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
     x2o3   =  2.0 / 3.0;

     satrec.init = 'y';
     satrec.t	 = 0.0;

     (
       satrec.no,
       method,
       ainv,  ao,    satrec.con41,  con42, cosio,
       cosio2,eccsq, omeosq, posq,
       rp,    rteosq,sinio , satrec.gsto,
       ) = _initl(
           satn, whichconst, satrec.ecco, epoch, satrec.inclo, satrec.no, satrec.method,
           satrec.afspc_mode
         );
     satrec.error = 0;

     """
     // sgp4fix remove this check as it is unnecessary
     // the mrt check in sgp4 handles decaying satellite cases even if the starting
     // condition is below the surface of te earth
//     if (rp < 1.0)
//       {
//         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
//         satrec.error = 5;
//       }
     """

     if omeosq >= 0.0 or satrec.no >= 0.0:

         satrec.isimp = 0;
         if rp < 220.0 / radiusearthkm + 1.0:
             satrec.isimp = 1;
         sfour  = ss;
         qzms24 = qzms2t;
         perige = (rp - 1.0) * radiusearthkm;

         #  - for perigees below 156 km, s and qoms2t are altered -
         if perige < 156.0:

             sfour = perige - 78.0;
             if perige < 98.0:
                 sfour = 20.0;
             #  sgp4fix use multiply for speed instead of pow
             qzms24temp =  (120.0 - sfour) / radiusearthkm;
             qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
             sfour  = sfour / radiusearthkm + 1.0;

         pinvsq = 1.0 / posq;

         tsi  = 1.0 / (ao - sfour);
         satrec.eta  = ao * satrec.ecco * tsi;
         etasq = satrec.eta * satrec.eta;
         eeta  = satrec.ecco * satrec.eta;
         psisq = fabs(1.0 - etasq);
         coef  = qzms24 * pow(tsi, 4.0);
         coef1 = coef / pow(psisq, 3.5);
         cc2   = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *
                        (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * satrec.con41 *
                        (8.0 + 3.0 * etasq * (8.0 + etasq)));
         satrec.cc1   = satrec.bstar * cc2;
         cc3   = 0.0;
         if satrec.ecco > 1.0e-4:
             cc3 = -2.0 * coef * tsi * j3oj2 * satrec.no * sinio / satrec.ecco;
         satrec.x1mth2 = 1.0 - cosio2;
         satrec.cc4    = 2.0* satrec.no * coef1 * ao * omeosq * \
                           (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *
                           (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
                           (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
                           (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
                           (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
         satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                        (etasq + eeta) + eeta * etasq);
         cosio4 = cosio2 * cosio2;
         temp1  = 1.5 * j2 * pinvsq * satrec.no;
         temp2  = 0.5 * temp1 * j2 * pinvsq;
         temp3  = -0.46875 * j4 * pinvsq * pinvsq * satrec.no;
         satrec.mdot     = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 * \
                            temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
         satrec.argpdot  = (-0.5 * temp1 * con42 + 0.0625 * temp2 *
                             (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                             temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4));
         xhdot1            = -temp1 * cosio;
         satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                              2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
         xpidot            =  satrec.argpdot+ satrec.nodedot;
         satrec.omgcof   = satrec.bstar * cc3 * cos(satrec.argpo);
         satrec.xmcof    = 0.0;
         if satrec.ecco > 1.0e-4:
             satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
         satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
         satrec.t2cof   = 1.5 * satrec.cc1;
         #  sgp4fix for divide by zero with xinco = 180 deg
         if fabs(cosio+1.0) > 1.5e-12:
             satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
         else:
             satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
         satrec.aycof   = -0.5 * j3oj2 * sinio;
         #  sgp4fix use multiply for speed instead of pow
         delmotemp = 1.0 + satrec.eta * cos(satrec.mo);
         satrec.delmo   = delmotemp * delmotemp * delmotemp;
         satrec.sinmao  = sin(satrec.mo);
         satrec.x7thm1  = 7.0 * cosio2 - 1.0;

         #----------- set variables if not deep space -----------
         if satrec.isimp != 1:

           cc1sq          = satrec.cc1 * satrec.cc1;
           satrec.d2    = 4.0 * ao * tsi * cc1sq;
           temp           = satrec.d2 * tsi * satrec.cc1 / 3.0;
           satrec.d3    = (17.0 * ao + sfour) * temp;
           satrec.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * \
                            satrec.cc1;
           satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
           satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
                            (12.0 * satrec.d2 + 10.0 * cc1sq));
           satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
                            12.0 * satrec.cc1 * satrec.d3 +
                            6.0 * satrec.d2 * satrec.d2 +
                            15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
     sgp4(satrec, 0.0);

     satrec.init = 'n';

     # sgp4fix return boolean. satrec.error contains any error codes
     return true;

@jit(cache=True)
def sgp4(satrec, tsince, whichconst=None):

     mrt = 0.0
     if whichconst is None:
          whichconst = satrec.whichconst
     
     temp4 =   1.5e-12;
     twopi = 2.0 * pi;
     x2o3  = 2.0 / 3.0;
     #  sgp4fix identify constants and allow alternate values
     tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
     vkmpersec     = radiusearthkm * xke/60.0;

     #  --------------------- clear sgp4 error flag -----------------
     satrec.t     = tsince;
     satrec.error = 0;
     satrec.error_message = None

     #  ------- update for secular gravity and atmospheric drag -----
     xmdf    = satrec.mo + satrec.mdot * satrec.t;
     argpdf  = satrec.argpo + satrec.argpdot * satrec.t;
     nodedf  = satrec.nodeo + satrec.nodedot * satrec.t;
     argpm   = argpdf;
     mm      = xmdf;
     t2      = satrec.t * satrec.t;
     nodem   = nodedf + satrec.nodecf * t2;
     tempa   = 1.0 - satrec.cc1 * satrec.t;
     tempe   = satrec.bstar * satrec.cc4 * satrec.t;
     templ   = satrec.t2cof * t2;

     if satrec.isimp != 1:

         delomg = satrec.omgcof * satrec.t;
         #  sgp4fix use mutliply for speed instead of pow
         delmtemp =  1.0 + satrec.eta * cos(xmdf);
         delm   = satrec.xmcof * \
                  (delmtemp * delmtemp * delmtemp -
                  satrec.delmo);
         temp   = delomg + delm;
         mm     = xmdf + temp;
         argpm  = argpdf - temp;
         t3     = t2 * satrec.t;
         t4     = t3 * satrec.t;
         tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 - \
                          satrec.d4 * t4;
         tempe  = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -
                          satrec.sinmao);
         templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +
                          satrec.t * satrec.t5cof);

     nm    = satrec.no;
     em    = satrec.ecco;
     inclm = satrec.inclo;
     if satrec.method == 'd':

         tc = satrec.t;
         (
             atime, em,    argpm,  inclm, xli,
             mm,    xni,   nodem,  dndt,  nm,
         ) = _dspace(
               satrec.irez,
               satrec.d2201, satrec.d2211, satrec.d3210,
               satrec.d3222, satrec.d4410, satrec.d4422,
               satrec.d5220, satrec.d5232, satrec.d5421,
               satrec.d5433, satrec.dedt,  satrec.del1,
               satrec.del2,  satrec.del3,  satrec.didt,
               satrec.dmdt,  satrec.dnodt, satrec.domdt,
               satrec.argpo, satrec.argpdot, satrec.t, tc,
               satrec.gsto, satrec.xfact, satrec.xlamo,
               satrec.no, satrec.atime,
               em, argpm, inclm, satrec.xli, mm, satrec.xni,
               nodem, nm
             );

     if nm <= 0.0:

         satrec.error_message = ('mean motion {0:f} is less than zero'
                                 .format(nm))
         satrec.error = 2;
         #  sgp4fix add return
         return false, false;

     am = pow((xke / nm),x2o3) * tempa * tempa;
     nm = xke / pow(am, 1.5);
     em = em - tempe;

     #  fix tolerance for error recognition
     #  sgp4fix am is fixed from the previous nm check
     if em >= 1.0 or em < -0.001:  # || (am < 0.95)

         satrec.error_message = ('mean eccentricity {0:f} not within'
                                 ' range 0.0 <= e < 1.0'.format(em))
         satrec.error = 1;
         #  sgp4fix to return if there is an error in eccentricity
         return false, false;

     #  sgp4fix fix tolerance to avoid a divide by zero
     if em < 1.0e-6:
         em  = 1.0e-6;
     mm     = mm + satrec.no * templ;
     xlm    = mm + argpm + nodem;
     emsq   = em * em;
     temp   = 1.0 - emsq;

     nodem  = nodem % twopi if nodem >= 0.0 else -(-nodem % twopi)
     argpm  = argpm % twopi
     xlm    = xlm % twopi
     mm     = (xlm - argpm - nodem) % twopi

     #  ----------------- compute extra mean quantities -------------
     sinim = sin(inclm);
     cosim = cos(inclm);

     #  -------------------- add lunar-solar periodics --------------
     ep     = em;
     xincp  = inclm;
     argpp  = argpm;
     nodep  = nodem;
     mp     = mm;
     sinip  = sinim;
     cosip  = cosim;
     if satrec.method == 'd':

         ep, xincp, nodep, argpp, mp = _dpper(
               satrec, satrec.inclo,
               'n', ep, xincp, nodep, argpp, mp, satrec.afspc_mode
             );
         if xincp < 0.0:

             xincp  = -xincp;
             nodep = nodep + pi;
             argpp  = argpp - pi;

         if ep < 0.0 or ep > 1.0:

             satrec.error_message = ('perturbed eccentricity {0:f} not within'
                                     ' range 0.0 <= e <= 1.0'.format(ep))
             satrec.error = 3;
             #  sgp4fix add return
             return false, false;

     #  -------------------- long period periodics ------------------
     if satrec.method == 'd':

         sinip =  sin(xincp);
         cosip =  cos(xincp);
         satrec.aycof = -0.5*j3oj2*sinip;
         #  sgp4fix for divide by zero for xincp = 180 deg
         if fabs(cosip+1.0) > 1.5e-12:
             satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
         else:
             satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;

     axnl = ep * cos(argpp);
     temp = 1.0 / (am * (1.0 - ep * ep));
     aynl = ep* sin(argpp) + temp * satrec.aycof;
     xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl;

     #  --------------------- solve kepler's equation ---------------
     u    = (xl - nodep) % twopi
     eo1  = u;
     tem5 = 9999.9;
     ktr = 1;
     #    sgp4fix for kepler iteration
     #    the following iteration needs better limits on corrections
     while fabs(tem5) >= 1.0e-12 and ktr <= 10:

         sineo1 = sin(eo1);
         coseo1 = cos(eo1);
         tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
         tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
         if fabs(tem5) >= 0.95:
             tem5 = 0.95 if tem5 > 0.0 else -0.95;
         eo1    = eo1 + tem5;
         ktr = ktr + 1;

     #  ------------- short period preliminary quantities -----------
     ecose = axnl*coseo1 + aynl*sineo1;
     esine = axnl*sineo1 - aynl*coseo1;
     el2   = axnl*axnl + aynl*aynl;
     pl    = am*(1.0-el2);
     if pl < 0.0:

         satrec.error_message = ('semilatus rectum {0:f} is less than zero'
                                 .format(pl))
         satrec.error = 4;
         #  sgp4fix add return
         return false, false;

     else:

         rl     = am * (1.0 - ecose);
         rdotl  = sqrt(am) * esine/rl;
         rvdotl = sqrt(pl) / rl;
         betal  = sqrt(1.0 - el2);
         temp   = esine / (1.0 + betal);
         sinu   = am / rl * (sineo1 - aynl - axnl * temp);
         cosu   = am / rl * (coseo1 - axnl + aynl * temp);
         su     = atan2(sinu, cosu);
         sin2u  = (cosu + cosu) * sinu;
         cos2u  = 1.0 - 2.0 * sinu * sinu;
         temp   = 1.0 / pl;
         temp1  = 0.5 * j2 * temp;
         temp2  = temp1 * temp;

         #  -------------- update for short period periodics ------------
         if satrec.method == 'd':

             cosisq                 = cosip * cosip;
             satrec.con41  = 3.0*cosisq - 1.0;
             satrec.x1mth2 = 1.0 - cosisq;
             satrec.x7thm1 = 7.0*cosisq - 1.0;

         mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + \
                 0.5 * temp1 * satrec.x1mth2 * cos2u;
         su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
         xnode = nodep + 1.5 * temp2 * cosip * sin2u;
         xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
         mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / xke;
         rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +
                 1.5 * satrec.con41) / xke;

         #  --------------------- orientation vectors -------------------
         sinsu =  sin(su);
         cossu =  cos(su);
         snod  =  sin(xnode);
         cnod  =  cos(xnode);
         sini  =  sin(xinc);
         cosi  =  cos(xinc);
         xmx   = -snod * cosi;
         xmy   =  cnod * cosi;
         ux    =  xmx * sinsu + cnod * cossu;
         uy    =  xmy * sinsu + snod * cossu;
         uz    =  sini * sinsu;
         vx    =  xmx * cossu - cnod * sinsu;
         vy    =  xmy * cossu - snod * sinsu;
         vz    =  sini * cossu;

         #  --------- position and velocity (in km and km/sec) ----------
         _mr = mrt * radiusearthkm
         r = (_mr * ux, _mr * uy, _mr * uz)
         v = ((mvt * ux + rvdot * vx) * vkmpersec,
              (mvt * uy + rvdot * vy) * vkmpersec,
              (mvt * uz + rvdot * vz) * vkmpersec)

     #  sgp4fix for decaying satellites
     if mrt < 1.0:

         satrec.error_message = ('mrt {0:f} is less than 1.0 indicating'
                                 ' the satellite has decayed'.format(mrt))
         satrec.error = 6;
         return false, false;

     return r, v;

def _gstime(jdut1):

     tut1 = (jdut1 - 2451545.0) / 36525.0;
     temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + \
             (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841;  #  sec
     temp = (temp * deg2rad / 240.0) % twopi # 360/86400 = 1/240, to deg, to rad

     #  ------------------------ check quadrants ---------------------
     if temp < 0.0:
         temp += twopi;

     return temp;

def getgravconst():
#  ------------ wgs-84 constants ------------
mu     = 398600.5;            #  in km3 / s2
radiusearthkm = 6378.137;     #  km
xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
tumin  = 1.0 / xke;
j2     =   0.00108262998905;
j3     =  -0.00000253215306;
j4     =  -0.00000161098761;
j3oj2  =  j3 / j2;

return tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2