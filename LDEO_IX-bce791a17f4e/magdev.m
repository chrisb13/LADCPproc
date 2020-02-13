function  [dev,d,h,i,f,x,y,z]=magdev(flat,flon,elevkm,nmax,igrf);
% function [dev,d,h,i,f,x,y,z]=magdev(flat,flon,elevkm,nmax,igrf);
% 
% compute magnetic deviation
% input:
%      flat    latitude (degree)
%      flon    longitude (degree)
%      elevkm  elevation above mean geoid (km)
%      nmax    number of harmonics
%      igrf    (need file, default 'IGRF95')
%
% output:
%     dev     mag. deviation (degree)
% 
% 
% based of FORTRAN ROUTINE GEOMAG.FOR
% more info under  http://fdd.gsfc.nasa.gov/IGRF.html
%
% M. Visbeck, LDEO FEB 2000
%

if nargin<5
 igrf='IGRF95';
 igrf='IGRF00';
end
eval(['gh=',igrf,';'])

if nargin<4
 nmax=10;
end

if nargin<3
 elevkm=0;
end



[x,y,z] = shval3(flat,flon,elevkm,nmax,gh);
[d,i,h,f]=dihf (x,y,z);

dev=d*180/pi;

if nargout<1
 disp([' model ',igrf,'  harmonics: ',int2str(nmax)])
 disp([' lat: ',num2str(flat)])
 disp([' lon: ',num2str(flon)])
 disp([' mag dev [deg]: ',num2str(d*180/pi)])
end

%=====================================================================

function [x,y,z] = shval3(flat,flon,elevkm,nmax,gh)
% function [x,y,z] = shval3(flat,flon,elevkm,nmax,gh)
% ================================================================
%
%       version 1.01
%
%       calculates field components from spherical harmonic (sh)
%       models.
%
%       input:
%           flat  - north latitude, in degrees
%           flon  - east longitude, in degrees
%           elevkm  - elevation above mean sea level 
%           nmax  - maximum degree and order of coefficients
%           gh    - schmidt quasi-normal internal spherical
%                   harmonic coefficients
%           iext  - external coefficients flag (= 0 if none)
%           ext   - the three 1st-degree external coefficients
%                   (not used if iext = 0)
%
%       output:
%           x     -  northward component
%           y     -  eastward component
%           z     -  vertically-downward component
%
%       based on subroutine 'igrf' by d. r. barraclough and
%       s. r. c. malin, report no. 71/1, institute of geological
%       sciences, u.k.
%
%       norman w. peddie, u.s. geological survey, mail stop 964,
%       federal center, box 25046, denver, colorado 80225
%
% ================================================================
%       the required sizes of the arrays used in this subroutine
%       depend on the value of nmax.  the minimum dimensions
%       needed are indicated in the table below.  (note that this
%       version is dimensioned for nmax of 14 or less).
%
% ================================================================
dtr = pi/180;
r=elevkm;
erad=6371.2;
a2=40680925.;
b2=40408588.;



slat = sin(flat*dtr);
aa = min(89.999,max(-89.999,flat));
clat = cos(aa*dtr);
sl(1) = sin(flon*dtr);
cl(1) = cos(flon*dtr);

x=0.0;
y=0.0;
z=0.0;
sd = 0.0;
cd = 1.0;
n=0;
l=1;
m=1;

npq = (nmax*(nmax+3))/2;
aa = a2*clat*clat;
bb = b2*slat*slat;
cc = aa+bb;
dd = sqrt(cc);
r=sqrt(elevkm*(elevkm+2.0*dd)+(a2*aa+b2*bb)/cc);
cd = (elevkm+dd)/r;
sd = (a2-b2)/dd*slat*clat/r;
aa = slat;
slat = slat*cd-clat*sd;
clat = clat*cd+aa*sd;
ratio = erad/r;
aa = sqrt(3.0);
p(1) = 2.0*slat;
p(2) = 2.0*clat;
p(3) = 4.5*slat*slat-1.5;
p(4) = 3.0*aa*clat*slat;
q(1) = -clat;
q(2) = slat;
q(3) = -3.0*clat*slat;
q(4) = aa*(slat*slat-clat*clat);


for k = 1: npq

if (n<m)
    m=0;
    n=n+1;
    rr = ratio^(n+2);
    fn=n;
end;

fm=m;

if (k>=5)
    if (m==n)
         aa = sqrt(1.0-.5/fm);
         j=k-n-1;
         p(k) = (1.0+1.0/fm)*aa*clat*p(j);
         q(k) = aa*(clat*q(j)+slat/fm*p(j));
         sl(m) = sl(m-1)*cl(1)+cl(m-1)*sl(1);
         cl(m) = cl(m-1)*cl(1)-sl(m-1)*sl(1);
    else 
         aa = sqrt(fn*fn-fm*fm);
         bb = sqrt((fn-1.0)^2-fm*fm)/aa;
         cc = (2.0*fn-1.0)/aa;
         i=k-n;
         j=k-2*n+1;
         p(k) = (fn+1.0)*(cc*slat/fn*p(i)-bb/(fn-1.0)*p(j));
         q(k) = cc*(slat*q(i)-clat/fn*p(i))-bb*q(j);
    end;
end;

aa = rr*gh(l);

if (m==0)
   x=x+aa*q(k);
   z=z-aa*p(k);
   l=l+1;
else 
   bb = rr*gh(l+1);
   cc = aa*cl(m)+bb*sl(m);
   x=x+cc*q(k);
   z=z-cc*p(k);
   if (clat>0.0)
     y=y+(aa*sl(m)-bb*cl(m))*fm*p(k)/((fn+1.0)*clat);
   else 
     y=y+(aa*sl(m)-bb*cl(m))*q(k)*slat;
   end;
   l=l+2;
end;

m=m+1;

end;

aa=x;
x=x*cd+z*sd;
z=z*cd-aa*sd;
          
return;

%==================================================================
function [d,i,h,f] = dihf(x,y,z)
% function [d,i,h,f] = dihf(x,y,z)
% ===============================================================
%
%       version 1.01
%
%       computes the geomagnetic elements d, i, h, and f from
%       x, y, and z.
%
%       input:
%           x   - northward component
%           y   - eastward component
%           z   - vertically-downward component
%
%       output:
%           d   - declination
%           i   - inclination
%           h   - horizontal intensity
%           f   - total intensity
%
%       a. zunde
%       usgs, ms 964, box 25046 federal center, denver, co  80225
%
% ===============================================================
sn=0.0001;

% ---------------------------------------------------------------
%       if d and i cannot be determined, set equal to NaN
% ---------------------------------------------------------------
h=sqrt(x*x+y*y);
f=sqrt(x*x+y*y+z*z);
if (f<sn)
  d=NaN;
  i=NaN;
else 
  i=atan2(z,h);
  if (h<sn)
    d=NaN;
  else 
    hpx = h+x;
    if (hpx<200)
      d=pi;
    else 
      d=2.0*atan2(y,hpx);
    end;
  end;
end;
return;
