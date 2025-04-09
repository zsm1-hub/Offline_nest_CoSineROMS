clear
fname='rst.seddy2.v8b19_2011_2013typhoon.nc';
ncdisp(fname,'/','full')
nc=netcdf(fname,'r')
Vstretching=nc{'Vstretching'}(:)

time=datestr(nc{'ocean_time'}(:)./86400+datenum(1900,1,1))
nc{'hc'}(:)
% 'alkalinity','TIC','oxygen','DDCA','CSDC','CLDC','SDOC','SDON','LDOC LDON',...
%     'DDSi','DD_C','DD_N','BAC_','Z2_C','Z2_N','Z1_C','Z1_N','S3CH','S3_C','S3_N',...
%     'S2CH','S2_C','S2_N','S1CH','S1_C','S1_N','PO4','SiOH4','NH4','NO3' 

fname='ini2012_L1.nc'
ncdisp(fname,'/','full')

nc=netcdf(fname,'r')
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
Vstretching=nc{'Vstretching'}(:)
NO3=squeeze(nc{'NO3'}(:,end,:,:));

NO3(abs(NO3)>1e20)=nan;
contourf(lon,lat,NO3);colorbar;