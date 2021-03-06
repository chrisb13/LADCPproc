function [gh,G,H] = IGRF95
% MATLAB routine to load Schmidt-normalized coefficients
% retrieved from ftp://nssdc.gsfc.nasa.gov/pub/models/igrf/
% igrf95.dat                  1 Kb    Mon Nov 13 00:00:00 1995
% ? C.E. Barton, Revision of International Geomagnetic Reference
% Field Released, EOS Transactions 77, #16, April 16, 1996.
% The coefficients are from the 1995 International Geomagnetic Reference Field
% Carlos Roithmayr, Jan. 22, 1997.
%++++++++++++++++++++++++++++++++++++++++++
% The number 1 is added to ALL subscripts since MATLAB can't have an array
% index of 0.  Units of Tesla
%
%  MARTIN VISBECK, LDEO Feb 2000
% !!! CAUTION  when updating the values NEVER us a 0. where use 0.1 instead!!!
% I use the zeros to throw unused coeffs away....!

G(2,1) = -29682e-9;
G(2,2) = -1789e-9; H(2,2) =  5318e-9;
G(3,1) = -2197e-9; H(3,1) =     0.0;
G(3,2) =  3074e-9; H(3,2) =   -2356e-9;
G(3,3) =  1685e-9; H(3,3) =  -425e-9;
G(4,1) =  1329e-9; H(4,1) =     0.0;
G(4,2) = -2268e-9; H(4,2) =  -263e-9;
G(4,3) =  1249e-9; H(4,3) =   302e-9;
G(4,4) =   769e-9; H(4,4) =  -406e-9;
G(5,1) =   941e-9; H(5,1) =      .0;
G(5,2) =   782e-9; H(5,2) =   262e-9;
G(5,3) =   291e-9; H(5,3) =  -232e-9;
G(5,4) =  -421e-9; H(5,4) =    98e-9;
G(5,5) =   116e-9; H(5,5) =  -301e-9;
G(6,1) =  -210e-9; H(6,1) =      .0;
G(6,2) =   352e-9; H(6,2) =    44e-9;
G(6,3) =   237e-9; H(6,3) =   157e-9;
G(6,4) =  -122e-9; H(6,4) =  -152e-9;
G(6,5) =  -167e-9; H(6,5) =   -64e-9;
G(6,6) =   -26e-9; H(6,6) =    99e-9;
G(7,1) =    66e-9; H(7,1) =      .0;
G(7,2) =    64e-9; H(7,2) =   -16e-9;
G(7,3) =    65e-9; H(7,3) =    77e-9;
G(7,4) =  -172e-9; H(7,4) =    67e-9;
G(7,5) =     2e-9; H(7,5) =   -57e-9;
G(7,6) =    17e-9; H(7,6) =     4e-9;
G(7,7) =   -94e-9; H(7,7) =    28e-9;
G(8,1) =    78e-9; H(8,1) =     -.0;
G(8,2) =   -67e-9; H(8,2) =   -77e-9;
G(8,3) =     1e-9; H(8,3) =   -25e-9;
G(8,4) =    29e-9; H(8,4) =     3e-9;
G(8,5) =     4e-9; H(8,5) =    22e-9;
G(8,6) =     8e-9; H(8,6) =    16e-9;
G(8,7) =    10e-9; H(8,7) =   -23e-9;
G(8,8) =    -2e-9; H(8,8) =    -3e-9;
G(9,1) =    24e-9; H(9,1) =      .0;
G(9,2) =     4e-9; H(9,2) =    12e-9;
G(9,3) =    -1e-9; H(9,3) =   -20e-9;
G(9,4) =    -9e-9; H(9,4) =     7e-9;
G(9,5) =   -14e-9; H(9,5) =   -21e-9;
G(9,6) =     4e-9; H(9,6) =    12e-9;
G(9,7) =     5e-9; H(9,7) =    10e-9;
G(9,8) =     0.1e-9; H(9,8) =   -17e-9;
G(9,9) =    -7e-9; H(9,9) =   -10e-9;
G(10,1) =     4e-9; H(10,1) =      .0;
G(10,2) =     9e-9; H(10,2) =   -19e-9;
G(10,3) =     1e-9; H(10,3) =    15e-9;
G(10,4) =   -12e-9; H(10,4) =    11e-9;
G(10,5) =     9e-9; H(10,5) =    -7e-9;
G(10,6) =    -4e-9; H(10,6) =    -7e-9;
G(10,7) =    -2e-9; H(10,7) =     9e-9;
G(10,8) =     7e-9; H(10,8) =     7e-9;
G(10,9) =     0.1e-9; H(10,9) =    -8e-9;
G(10,10) =    -6e-9; H(10,10) =     1e-9;
G(11,1) =    -3e-9; H(11,1) =      .0;
G(11,2) =      -4e-9; H(11,2) =     2e-9;
G(11,3) =       2e-9; H(11,3) =     1e-9;
G(11,4) =      -5e-9; H(11,4) =     3e-9;
G(11,5) =      -2e-9; H(11,5) =     6e-9;
G(11,6) =       4e-9; H(11,6) =    -4e-9;
G(11,7) =       3e-9; H(11,7) =     0.1e-9;
G(11,8) =       1e-9; H(11,8) =    -2e-9;
G(11,9) =       3e-9; H(11,9) =     3e-9;
G(11,10) =      3e-9; H(11,10) =    -1e-9;
G(11,11) =      0.1e-9; H(11,11) =    -6e-9;

% prepare compressed array

g=reshape(G',11*11,1)*1e9;
h=reshape(H',11*11,1)*1e9;

gh=reshape([g,h]',2*11*11,1);
% here is where zeros have a meaning....!!!!
ii=find(gh==0);
gh(ii)=[];

return;


