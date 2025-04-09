%
%  D_ROMS2ROMS:  Driver script to create a ROMS initial conditions
%
%  This a user modifiable script that can be used to prepare ROMS
%  initial conditions NetCDF file from another ROMS application. It
%  sets-up all the necessary parameters and variables. USERS can use
%  this as a prototype for their application.
%

% svn $Id: d_roms2roms.m 1156 2023-02-18 01:44:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%==========================================================================
% This script creates initial conditions from ROMS NWA dataset.
%==========================================================================

% Set file names.

clear all;close all;clc
addpath('D:\LIN2023\model\Oforc_OGCM\DATA\WOAPISCES')
addpath('D:\LIN2023\model\Oforc_OGCM\DATA\WOA2009')
addpath('../nest/Preprocessing_tools/')
addpath(genpath('D:\LIN2023\crocotools\'))
addpath(genpath('D:\UTILITIES/'))
addpath(genpath('D:\LIN2023\model\offline'))
addpath(genpath('F:\offline\L1\'))
addpath('D:\LIN2023\model\')
addpath('F:\NSCS')
% '01-Jan-2012 06:00:00 tindex=2

NWAdata = 'F:\offline\L1\rst.seddy2.s8b01_L1.nc';
NWAgrid = 'roms_grd_modi_sbry.nc';
datestr(ncread(NWAdata,'ocean_time')./86400+datenum(1900,1,1))

GRDname = 'roms_grd_modi_sbry_0921.nc.1';
INIname = 'F:\offline\L2\ini2012_L2_0702_01.nc';

CREATE = true;                   % logical switch to create NetCDF
report = false;                  % report vertical grid information

IniRec = 4;                      % NWA time record for initialization

time=datestr(ncread(NWAdata,'ocean_time')./86400+datenum(1900,1,1))
%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

[Lr,Mr] = size(nc_read(GRDname,'h'));

Lu = Lr-1;
Lv = Lr;
Mu = Mr;
Mv = Mr-1;

S.ncname      = INIname;    % output NetCDF file

S.spherical   = 1;          % spherical grid

S.Lm          = Lr-2;       % number of interior RHO-points, X-direction
S.Mm          = Mr-2;       % number of interior RHO-points, Y-direction
S.N           = 30;         % number of vertical levels at RHO-points
S.NT          = 2;          % total number of tracers

S.Vtransform  = 2;          % vertical transfomation equation
S.Vstretching = 4;          % vertical stretching function

S.theta_s     = 7.0;        % S-coordinate surface control parameter
S.theta_b     = 2.0;        % S-coordinate bottom control parameter
S.Tcline      = 500;      % S-coordinate surface/bottom stretching width
S.hc          = S.Tcline;   % S-coordinate stretching width

%--------------------------------------------------------------------------
% Set variables to process.
%--------------------------------------------------------------------------

VarGrd  = {'spherical',                                               ...
           'Vtransform', 'Vstretching',                               ...
           'theta_s', 'theta_b', 'Tcline', 'hc',                      ...
           's_rho', 'Cs_r', 's_w', 'Cs_w', 'h'};

if (S.spherical),
  VarGrd = [VarGrd, 'lon_rho', 'lat_rho',                             ...
                    'lon_u', 'lat_u', 'lon_v', 'lat_v'];
else
  VarGrd = [VarGrd, 'x_rho', 'y_rho',                                 ...
                    'x_u', 'y_u', 'x_v', 'y_v'];
end

VarIni = {'zeta', 'ubar', 'vbar', 'u', 'v', 'temp', 'salt',...
    'alkalinity','TIC','oxygen','DDCA','CSDC','CLDC','SDOC','SDON','LDOC','LDON',...
    'DDSi','DD_C','DD_N','BAC_','Z2_C','Z2_N','Z1_C','Z1_N','S3CH','S3_C','S3_N',...
    'S2CH','S2_C','S2_N','S1CH','S1_C','S1_N','PO4','SiOH4','NH4','NO3' };

%  Set intepolation parameters.

method = 'linear';             % linear interpolation
offset = 10;                   % number of extra points for sampling
RemoveNaN = true;              % remove NaN with nearest-neighbor
Rvector = true;                % interpolate vectors to RHO-points

%--------------------------------------------------------------------------
%  Get parent and target grids structures. The depths are for an
%  unperturbed state (zeta = 0).
%--------------------------------------------------------------------------

%  Get Parent grid structure, P.

P = get_roms_grid(NWAgrid, NWAdata);

%  Set surface-depths to zero to bound surface interpolation. This is
%  specific for this application.

N = P.N;

% P.z_r(:,:,N) = 0;
% P.z_u(:,:,N) = 0;
% P.z_v(:,:,N) = 0;

%  Get Target grid structure, T.

T = get_roms_grid(GRDname, S);

%  If vector rotation is required in the parent grid, interpolate
%  rotation angle (parent to target) and add it to target grid
%  structure.

T.parent_angle = roms2roms(NWAgrid, P, T, 'angle', [], Rvector,       ...
                           method, offset, RemoveNaN);

%--------------------------------------------------------------------------
%  Interpolate initial conditions from source data to application grid.
%--------------------------------------------------------------------------

disp(' ')
disp(['***********************************************************']);
disp(['** Interpolating initial conditions from NWA to GOM grid **']);
disp(['***********************************************************']);

%  The NWA data has a time coordinate in seconds, which starts
%  on 1-Jan-1900.

time = nc_read(NWAdata,'ocean_time',IniRec);
epoch = datenum('1-Jan-1900');
mydate = datestr(epoch+time/86400);

disp(' ')
disp(['** Processing: ',mydate,' **']);
disp(' ')

%  Set initial conditions time (seconds). The time coordinate for this
%  ROMS application is "seconds since 2000-01-01 00:00:00".

I.ocean_time = time;

%  Interpolate initial conditions.

for var = VarIni(1:3)
  field = char(var);
  I.(field) = roms2roms(NWAdata, P, T, field, IniRec, Rvector,        ...
                        method, offset, RemoveNaN);
end

%%%change by zsm ------define vertical coordinate by yourself!
Pzeta=ncread(NWAdata,'zeta');
Pzeta=squeeze(Pzeta(:,:,IniRec));
Pzr=zlevs([P.h]',[Pzeta]',P.theta_s,P.theta_b, ...
    P.hc,N,'r',P.Vtransform);
Pzw=zlevs([P.h]',[Pzeta]',P.theta_s,P.theta_b, ...
    P.hc,N,'w',P.Vtransform);
P.z_w = permute(Pzw,[3 2 1]);
P.z_r = permute(Pzr,[3 2 1]);
P.z_u = permute(rho2u_3d(Pzr),[3 2 1]);
P.z_v = permute(rho2v_3d(Pzr),[3 2 1]);

Tzr=zlevs([T.h]',[I.zeta]',T.theta_s,T.theta_b, ...
    T.hc,N,'r',T.Vtransform);
Tzw=zlevs([T.h]',[I.zeta]',T.theta_s,T.theta_b, ...
    T.hc,N,'w',T.Vtransform);
T.z_w = permute(Tzw,[3 2 1]);
T.z_r = permute(Tzr,[3 2 1]);
T.z_u = permute(rho2u_3d(Tzr),[3 2 1]);
T.z_v = permute(rho2v_3d(Tzr),[3 2 1]);
T.Hz = abs(T.z_w(:,:,2:end)-T.z_w(:,:,1:end-1));
clear Pzeta;clear Pzr;clear Tzr;clear Pzw;clear Tzw;

kgrid = 0;
[T.s_rho,T.Cs_r] = stretching(T.Vstretching,             ...
    T.theta_s,                 ...
    T.theta_b,                 ...
    T.hc, N,              ...
    kgrid, false);

kgrid = 1;
[T.s_w,  T.Cs_w] = stretching(T.Vstretching,             ...
    T.theta_s,                 ...
    T.theta_b,                 ...
   T.hc, N, kgrid);


%%%%%%%%%%%% end of zsm

for var = VarIni(4:end)
  field = char(var);
  I.(field) = roms2roms(NWAdata, P, T, field, IniRec, Rvector,        ...
                        method, offset, RemoveNaN);
end

%  Rotate interpolated 3D velocity at RHO-points to TRUE North and East.
%  Need to interpolate Parent grid rotation angle to Target grid.

irotate = 0;               % rotate for (XI,ETA) to (lon,lat)

[Urho,Vrho] = rotate_vec(I.u, I.v, T.parent_angle, irotate);

%  Rotate resulting 3D velocity (RHO-points) to target grid angle and
%  average to staggered C-grid locations.

[I.u,I.v] = roms_vectors(Urho, Vrho, T.angle, T.mask_u, T.mask_v);

%  Compute barotropic velocities by vertically integrating (u,v).

[I.ubar,I.vbar] = uv_barotropic(I.u, I.v, T.Hz);

%--------------------------------------------------------------------------
%  Create initial condition Netcdf file.
%--------------------------------------------------------------------------

if (CREATE),
  [status] = c_initial(S);

%  Set attributes for "ocean_time".

  avalue = 'seconds since 1900-01-01 00:00:00';
  [status] = nc_attadd(INIname,'units',avalue,'ocean_time');
  
  avalue = 'gregorian';
  [status] = nc_attadd(INIname,'calendar',avalue,'ocean_time');

%  Set global attributes.

  avalue = 'TWS, ~0.5 km resolution, Grid a';
  [status] = nc_attadd(INIname,'title',avalue);

  avalue = 'ROMS NWA application';
  [status] = nc_attadd(INIname,'data_source',avalue);

  [status] = nc_attadd(INIname,'grd_file',GRDname);
end,

%--------------------------------------------------------------------------
%  Write out initial conditions.
%--------------------------------------------------------------------------

if (CREATE),
  disp(' ')
  disp(['** Writing initial conditions **']);
  disp(' ')

  for var = VarGrd,
    field = char(var);
    [status] = nc_write(INIname, field, T.(field));
  end
  
  IniRec = 1;

  field = 'ocean_time';
  [err.(field)] = nc_write(INIname, field, I.(field), IniRec);
  %%%%%%%%%%%%%%%%% add by zsm --filled the 1e37 to the zeros
  I.zeta(I.zeta==0)=1e37;
  I.ubar(I.ubar==0)=1e37;
  I.vbar(I.vbar==0)=1e37;
  I.u(I.u==0)=1e37;
  I.v(I.v==0)=1e37;
  I.temp(I.temp==0)=1e37;
  I.salt(I.salt==0)=1e37;
  I.alkalinity(I.alkalinity==0)=1e37;
  I.TIC(I.TIC==0)=1e37;
  I.oxygen(I.oxygen==0)=1e37;
  I.DDCA(I.DDCA==0)=1e37;
  I.CSDC(I.CSDC==0)=1e37;
  I.CLDC(I.CLDC==0)=1e37;
  I.SDOC(I.SDOC==0)=1e37;
  I.SDON(I.SDON==0)=1e37;
  I.LDOC(I.LDOC==0)=1e37;
  I.LDON(I.LDON==0)=1e37;
  I.DDSi(I.DDSi==0)=1e37;
  I.DD_C(I.DD_C==0)=1e37;
  I.DD_N(I.DD_N==0)=1e37;
  I.BAC_(I.BAC_==0)=1e37;
  I.Z2_C(I.Z2_C==0)=1e37;
  I.Z2_N(I.Z2_N==0)=1e37;
  I.Z1_C(I.Z1_C==0)=1e37;
  I.Z1_N(I.Z1_N==0)=1e37;
  I.S3CH(I.S3CH==0)=1e37;
  I.S3_C(I.S3_C==0)=1e37;
  I.S3_N(I.S3_N==0)=1e37;
  I.S2CH(I.S2CH==0)=1e37;
  I.S2_C(I.S2_C==0)=1e37;
  I.S2_N(I.S2_N==0)=1e37;
  I.S1CH(I.S1CH==0)=1e37;
  I.S1_C(I.S1_C==0)=1e37;
  I.S1_N(I.S1_N==0)=1e37;
  I.PO4(I.PO4==0)=1e37;
  I.SiOH4(I.SiOH4==0)=1e37;
  I.NH4(I.NH4==0)=1e37;
  I.NO3(I.NO3==0)=1e37;
  % 'alkalinity','TIC','oxygen','DDCA','CSDC','CLDC','SDOC','SDON','LDOC','LDON',...
  %   'DDSi','DD_C','DD_N','BAC_','Z2_C','Z2_N','Z1_C','Z1_N','S3CH','S3_C','S3_N',...
  %   'S2CH','S2_C','S2_N','S1CH','S1_C','S1_N','PO4','SiOH4','NH4','NO3' 

  for var = VarIni
    field = char(var);
    [err.(field)] = nc_write(INIname, field, I.(field), IniRec);
  end
end