
clear all;close all;clc
% addpath('D:\LIN2023\model/nest/Preprocessing_tools/')
addpath(genpath('D:\LIN2023\crocotools\'))
addpath(genpath('D:\UTILITIES/'))
% addpath(genpath('D:\LIN2023\model\offline'))
addpath(genpath('D:\LIN2023\model\luwf\rutgers'))
addpath('F:\offline\L1')
addpath('D:\LIN2023\model\')
addpath('D:\LIN2023\model\luwf')

% diary off;

NWAdata = ['his.seddy2.s7b2_L1_10days_0039.nc'];
temp=ncread(NWAdata,'temp');
temp(abs(temp)>10000)=nan;
contourf(squeeze(temp(:,:,end,end)))
colorbar




clear all;close all;clc
% addpath('D:\LIN2023\model/nest/Preprocessing_tools/')
addpath(genpath('D:\LIN2023\crocotools\'))
addpath(genpath('D:\UTILITIES/'))
% addpath(genpath('D:\LIN2023\model\offline'))
addpath(genpath('D:\LIN2023\model\luwf\rutgers'))
addpath('F:\offline\L1')
addpath('D:\LIN2023\model\')
addpath('D:\LIN2023\model\luwf')

for files=6:12

    % NWAdata = ['F:\offline\L1\his.seddy2.s7b2_L1_10days_00',num2str(files+32,'%02d'),'.nc'];
    NWAdata = ['F:\offline\L1\avg.seddy2.s8b01_L1_00',num2str(files,'%02d'),'.nc'];
    fname=['D:\LIN2023\model\luwf\offlinetools\bry_s7b2_L2_',num2str(files),'.nc']
    % ncdisp(fname)
    % ncread(fname,'bry_time')
    time=ncread(NWAdata,'ocean_time');
    % T=time(:)-(datenum(2012,1,1)-datenum(1900,1,1)).*86400;
    T=(time(:)./86400+datenum(1900,1,1)).*86400;
    datestr(T./86400)

    nw=netcdf(fname,'w')
    ncdisp(fname)
    % nw{'bry_time'}(:)=time-86400.*(datenum(2012,6,09)-datenum(1900,1,1))-360;
    % nw{'bry_time'}(:)=time-86400.*(datenum(2012,7,1)-datenum(1900,1,1))-22*3600-12*60;
    nw{'bry_time'}(:)=time-86400.*(datenum(2012,7,1)-datenum(1900,1,1))-24*3600-3600;
    close(nw)

    nc = netcdf(fname,'write');
    result = redef(nc);
    s_rho=nc{'s_rho'}(:);
    % eta_rho=442;
    % xi_rho=597;
    bry_time=nc{'bry_time'}(:);
    nc('xi_rho')=597;
    nc('eta_rho')=442;
    %
    %  Create variables
    %
    nc{'NO_3_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'NO_3_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'NO_3_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'NO_3_west'} = ncdouble('bry_time','s_rho','eta_rho') ;


    nc{'PO_4_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'PO_4_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'PO_4_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'PO_4_west'} = ncdouble('bry_time','s_rho','eta_rho') ;


    nc{'TInC_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'TInC_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'TInC_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'TInC_west'} = ncdouble('bry_time','s_rho','eta_rho') ;



    nc{'TALK_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'TALK_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'TALK_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'TALK_west'} = ncdouble('bry_time','s_rho','eta_rho') ;


    % Oxyg

    nc{'Oxyg_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'Oxyg_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'Oxyg_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'Oxyg_west'} = ncdouble('bry_time','s_rho','eta_rho') ;


    % SiOH
    nc{'SiOH_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'SiOH_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'SiOH_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'SiOH_west'} = ncdouble('bry_time','s_rho','eta_rho') ;



    nc{'NH_4_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'NH_4_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'NH_4_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'NH_4_west'} = ncdouble('bry_time','s_rho','eta_rho') ;




    % BACT
    nc{'BACT_south'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'BACT_east'} = ncdouble('bry_time','s_rho','eta_rho') ;
    nc{'BACT_north'} = ncdouble('bry_time','s_rho','xi_rho') ;
    nc{'BACT_west'} = ncdouble('bry_time','s_rho','eta_rho') ;

    close(nc);

    ncdisp(fname)
    nc=netcdf(fname,'w')

    nc{'NO_3_south'}(:) =nc{'NO3_south'}(:);
    nc{'NO_3_east'}(:) =nc{'NO3_east'}(:);
    nc{'NO_3_north'}(:) =nc{'NO3_north'}(:);
    nc{'NO_3_west'}(:) =nc{'NO3_west'}(:);
    nc{'PO_4_south'}(:) =nc{'PO4_south'}(:);
    nc{'PO_4_east'}(:) =nc{'PO4_east'}(:);
    nc{'PO_4_north'}(:) =nc{'PO4_north'}(:);
    nc{'PO_4_west'}(:) =nc{'PO4_west'}(:);
    nc{'TInC_south'}(:) =nc{'TIC_south'}(:);
    nc{'TInC_east'}(:) =nc{'TIC_east'}(:);
    nc{'TInC_north'}(:) =nc{'TIC_north'}(:);
    nc{'TInC_west'}(:) =nc{'TIC_west'}(:);
    nc{'TALK_south'}(:) =nc{'alkalinity_south'}(:);
    nc{'TALK_east'}(:) =nc{'alkalinity_east'}(:);
    nc{'TALK_north'}(:) =nc{'alkalinity_north'}(:);
    nc{'TALK_west'}(:) =nc{'alkalinity_west'}(:);

    nc{'Oxyg_south'}(:) =nc{'oxygen_south'}(:);
    nc{'Oxyg_east'}(:) =nc{'oxygen_east'}(:);
    nc{'Oxyg_north'}(:)=nc{'oxygen_north'}(:);
    nc{'Oxyg_west'}(:) =nc{'oxygen_west'}(:);
    nc{'SiOH_south'}(:) =nc{'SiOH4_south'}(:);
    nc{'SiOH_east'}(:) =nc{'SiOH4_east'}(:);
    nc{'SiOH_north'}(:) =nc{'SiOH4_north'}(:);
    nc{'SiOH_west'}(:) =nc{'SiOH4_west'}(:);

    nc{'NH_4_south'}(:) =nc{'NH4_south'}(:);
    nc{'NH_4_east'}(:) =nc{'NH4_east'}(:);
    nc{'NH_4_north'}(:) =nc{'NH4_north'}(:);
    nc{'NH_4_west'}(:) =nc{'NH4_west'}(:);
    nc{'BACT_south'}(:) =nc{'BAC__south'}(:);
    nc{'BACT_east'}(:) =nc{'BAC__east'}(:);
    nc{'BACT_north'}(:) =nc{'BAC__north'}(:);
    nc{'BACT_west'}(:) =nc{'BAC__west'}(:);
    close(nc);

end

files=1
fname=['bry_s7b2_L2_',num2str(files),'.nc']
datestr(ncread(fname,'bry_time')./86400+datenum(1900,1,1))

aa=ncread(fname,'ubar_east');
bb=squeeze(aa(:,end,:));
% sed -i 's/DD_C_time/bry_time/g; s/BACT_time/bry_time/g; s/DD_N_time/bry_time/g' varinfo_bry_L2.dat
% 
% sed -i 's/DDCA_time/bry_time/g; s/DDSi_time/bry_time/g' varinfo_bry_L2.dat
% 
% sed -i 's/LDON_time/bry_time/g; s/LDOC_time/bry_time/g' varinfo_bry_L2.dat
% 
% sed -i 's/SDON_time/bry_time/g; s/SDOC_time/bry_time/g' varinfo_bry_L2.dat
% 
% sed -i 's/CLDC_time/bry_time/g; s/CSDC_time/bry_time/g' varinfo_bry_L2.dat
% 
% sed -i 's/Z1_C_time/bry_time/g; s/Z1_N_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/Z2_C_time/bry_time/g; s/Z2_N_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/S3_C_time/bry_time/g; s/S3_N_time/bry_time/g; s/S3CH_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/S2_C_time/bry_time/g; s/S2_N_time/bry_time/g; s/S2CH_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/S1_C_time/bry_time/g; s/S1_N_time/bry_time/g; s/S1CH_time/bry_time/g' varinfo_bry_L2.dat
% 
% 
% sed -i 's/SiOH_time/bry_time/g; s/NH_4_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/PO_4_time/bry_time/g; s/NO_3_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/TALK_time/bry_time/g; s/Oxyg_time/bry_time/g' varinfo_bry_L2.dat
% sed -i 's/TInC_time/bry_time/g' varinfo_bry_L2.dat

% 
% SiOH_time
% NH_4_time
% PO_4_time
% NO_3_time
% TALK_time
% Oxyg_time
% TInC_time