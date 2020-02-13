%CB: 14/01/2020

%icebergs ctds are:
%[5, 8, 9, 12, 13, 14, 15, 16, 17, 20, 23, 24, 27, 28, 31, 34, 35, 36, 37, 38, 41, 44, 46, 47, 49, 50, 51, 53, 55, 57, 58, 60, 61, 63]

% these ones process fine:
% for i = [5,8,9,12,14,20, 24, 28, 31, 34, 35, 36, 37, 44, 46, 49]

% problem processing 13, 15, 16, 17, 23, 27, 38, 41, 42

% Specifically: 13, 15, 16, 23, 41, 42, 50, 60, 61, 63
% Index exceeds the number of array elements (0).

% Specifically: 17, 27, 38, 51, 57, 58
% Non-finite values in d.izm --- try processing with p.getdepth == 1
% FIXED! using p.getdepth == 1

% 53 and 55 was different: 
% Error using prepinv (line 680)
% not enough data to process station

%%%%%%

% 47 seth problem
% 53 and 55 'not enough data to process station'
% 61 Cannot determine time offset between CTD and LADCP time series --- aborting

for i = [ 63]
       
    disp(['==> Processing_cast: ',int2str(i)])
    ofol=sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\processed\\JR19002_%03d',i);
    disp(['==> Creating outputfolder ',ofol]);
    system(append('mkdir ',ofol));
    ctdf=sprintf('%03d',i);
    disp(['==> Converting CTD to ascii file: ',ctdf])
    convert_CTD_to_ascii(i);
    process_cast(i);     
    ofiles=sprintf('C:\\Users\\Chris\\VBoxshared\\Chris_LADCP\\processing\\processed\\%03d*',i);
    movefile(ofiles,ofol);
    
    % create output file that can be read in python

    matout=sprintf('%03d',i);
    ofile_mat=append(ofol,'\',matout,'.mat');
    efile=append(ofol,'\',matout,'_pyexport.mat');
    
    ladcp=load(ofile_mat);
    u=ladcp.dr.u;
    v=ladcp.dr.v;
    z=ladcp.dr.z;
    press=ladcp.dr.p;
    lat=ladcp.dr.lat;
    lon=ladcp.dr.lon;
    date=ladcp.dr.date;
    ctdt=ladcp.dr.ctd_t;
    ctds=ladcp.dr.ctd_s;
    save(efile,'u','v','z','press','lat','lon','date','ctdt','ctds');
    
    disp(['==> Finished, processing_cast: ',int2str(i)])
    close all
    clear all
    
    
end