%% convert_CTD_to_ascii
%%
%% Takes 1Hz .mat files from CTD cast with NMEA GPS stream and converts to ASCII format for reading into LDEO IX
%% JR15006
%% Andrew Meijers
%% 8/4/16 - on the Antarctic circle
%% Ed AM 01/03/19 - just north of Ant circle, for JR18005

%% CB used on JR19002

function []=convert_CTD_to_ascii(stn)

if ~exist([sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\ctd\\JR19002_ctd_%03d.1Hz',stn)],'file')
    return 
else
end

ctdin=load([sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\ctd\\JR19002_ctd_%03d.1Hz',stn)],'-mat');

tm=[ctdin.gtime(1); ctdin.time_elapsed];
p=[ctdin.gtime(2);ctdin.press];
s=[ctdin.gtime(3);ctdin.salin];
t=[ctdin.gtime(4);ctdin.temp];
lat=[ctdin.gtime(5);ctdin.latscan];
lon=[ctdin.gtime(6);ctdin.lonscan];

CTDall=[tm,p,s,t,lat,lon];

flnm=sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\ctd\\JR19002_ctd_%03d.1Hz_ascii',stn);

save(flnm,'-ascii','-double','CTDall')


    