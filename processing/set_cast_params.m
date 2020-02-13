% 14/01/2020
% adapted from: p.cruise_id = 'JR18005'; p.whoami = 'A. Meijers & C. Pimm';

close all
more off;
% mkSADCP('N:\adcp','N:\adcp\SADCP.mat');

% the variable stn contains the station number (1st parameter to process cast.m)
% stn='string';

% The script set cast params.m is called (twice!) from process cast.m. 
% The following Matlab structures are typically modified in set cast params.m:
% f contains file- and path names; 
% p contains parameters used for reading and preparing the data files; 
% ps contains parameters used for calculating the absolute velocities.


f.res = sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\processed\\%03d',stn);
f.checkpoints = sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\checkpoints\\%03d',stn);
f.ladcpdo = sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\raw\\JR19002M_%03d.000',stn);
f.ladcpup = sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\raw\\JR19002S_%03d.000',stn);

p.saveplot = [1:4 11 13:14];
p.saveplot_pdf = [1:4 11 13:14];
% f.ctd

p.cruise_id = 'JR19002';
p.whoami = 'C. Bull';
p.ladcp_station = stn;

p.name = sprintf('%s case #%d',p.cruise_id,p.ladcp_station);

% To avoid this error:
% Non-finite values in d.izm --- try processing with p.getdepth == 1
if stn == 013
    p.getdepth = 1;
end

if stn == 017
    p.getdepth = 1;
end

if stn == 027
    p.getdepth = 1;
end

if stn == 038
    p.getdepth = 1;
end

if stn == 042
    p.getdepth = 1;
end

if stn == 050
    p.getdepth = 1;
end

if stn == 051
    p.getdepth = 1;
end

if stn == 053
    p.getdepth = 1;
end

% still ultimately crashes with: not enough data to process station
if stn == 055
    p.getdepth = 1;
end

if stn == 057
    p.getdepth = 1;
end

if stn == 058
    p.getdepth = 1;
end

if stn == 063
    p.getdepth = 1;
end

p.edit_mask_dn_bins = [];
p.edit_mask_up_bins = [];

p.checkpoints = [1];

f.ctd = sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\ctd\\JR19002_ctd_%03d.1Hz_ascii',stn);
% from andrew
f.ctd_header_lines = 1;
f.ctd_fields_per_line = 6;
f.ctd_pressure_field = 2;
f.ctd_temperature_field = 4;
f.ctd_salinity_field = 3;
f.ctd_time_field = 1;
f.ctd_time_base = 0;

f.nav = f.ctd;
f.nav_header_lines = 1;
f.nav_fields_per_line = 6;
f.nav_lat_field = 5;
f.nav_time_field = f.ctd_time_field;
f.nav_time_base = f.ctd_time_base;
f.nav_lon_field = 6;





% sadcp_dir='../SADCP/Matlab/Processing';

sadcp_dir='C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\sadcp\\Processed\\';
sadcp_files=dir(fullfile(sadcp_dir,'_JR18005_000_000000_*_abs.mat'));
if isempty(sadcp_files)
        warning('No SADCP data');
else
    f.sadcp=sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\sadcp\\JR18005_%03d_sadcp.mat',stn);
    % load fullfile(../SADCP/Matlab/Processing/JR310000_000000_1_abs.mat
    [~,adcp_ind]=max([sadcp_files.datenum]);
    fprintf(1,'Loading SADCP data from %s\n',...
        fullfile(sadcp_dir,sadcp_files(adcp_ind).name));
    load(fullfile(sadcp_dir,sadcp_files(adcp_ind).name));
    tim_sadcp=julian(datevec(OS75_abs.nav.txy1(1,:)+datenum(2019,1,1)));

    good_ind=find(~isnan(tim_sadcp));

    tim_sadcp=tim_sadcp(good_ind);
    lat_sadcp=OS75_abs.nav.txy1(3,good_ind)';
    lon_sadcp=OS75_abs.nav.txy1(2,good_ind)';
    u_sadcp=squeeze(OS75_abs.vel(:,1,good_ind));
    v_sadcp=squeeze(OS75_abs.vel(:,2,good_ind));
    z_sadcp=bindepth(:,good_ind(1));

    %%%%%%%
    %Mask any data that is below 25% good
    pg = [OS75_ave_ping.pg];
    pg2=pg(:,good_ind);
    mask_pg=ones(size(pg2));
    mask_pg(pg2<25) = NaN;

    u_sadcp = u_sadcp.*mask_pg;
    v_sadcp = v_sadcp.*mask_pg;
    %%%%%%%

    save (f.sadcp,'*_sadcp');
    clear *_sadcp OS75* bindepth adcp_ind
end

p.checkpoints = [];
