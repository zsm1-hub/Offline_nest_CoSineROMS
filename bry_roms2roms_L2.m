%
%  D_OBC_ROMS2ROMS:  Driver script to create a ROMS boundary conditions
%
%  This a user modifiable script that can be used to prepare ROMS open
%  boundary conditions NetCDF file from another ROMS application. It
%  sets-up all the necessary parameters and variables. USERS can use
%  this as a prototype for their application.
%

% svn $Id: d_obc_roms2roms.m 1156 2023-02-18 01:44:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%==========================================================================
% This script creates boundary conditions from ROMS NWA dataset.
%==========================================================================

% Set file names (Input data can be in an OpenDAP server).

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
NWAgrid = 'roms_grd_modi_sbry.nc';
datestr(ncread(NWAdata,'ocean_time')./86400+datenum(1900,1,1))

GRDname = 'roms_grd_modi_sbry_0921.nc.1';
BRYname = 'bry_L2_7.nc';

StrDay = datenum('9-Jul-2012 01:00');        % starting day to process
EndDay = datenum('14-Jul-2012 00:00');        % ending   day to process
REFTIME=datenum(2012,1,1)-datenum(1900,1,1);

CREATE = true;                          % logical switch to create NetCDF
WRITE  = true;                          % logical switch to write out data
report = false;                         % report vertical grid information

%  Set logical switch to compute time dependent vertical thichnesses (Hz)
%  when computing vertically integrated 2D momentum (ubar,vbar).

TimeDependent = true;

%  Initialize unlimited dimension record counter. Notice that if we want
%  to restart the computations, we can set CREATE = false and get the
% record of the last boundary conditions processed for appending.

if (CREATE)
  BryRec = 0;
else
  BryRec = length(nc_read(BRYname,'bry_time'));
end

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

[Lr,Mr] = size(nc_read(GRDname,'h'));

Lu = Lr-1;
Lv = Lr;
Mu = Mr;
Mv = Mr-1;

S.ncname      = BRYname;    % output NetCDF file

S.spherical   = 1;          % spherical grid

S.Lm          = Lr-2;       % number of interior RHO-points, X-direction
S.Mm          = Mr-2;       % number of interior RHO-points, Y-direction
S.N           = 30;         % number of vertical levels at RHO-points
S.NT          = 33;          % total number of tracers

S.Vtransform  = 2;          % vertical transfomation equation
S.Vstretching = 4;          % vertical stretching function

S.theta_s     = 7.0;        % S-coordinate surface control parameter
S.theta_b     = 2.0;        % S-coordinate bottom control parameter
S.Tcline      = 500.0;      % S-coordinate surface/bottom stretching width
S.hc          = S.Tcline;   % S-coordinate stretching width

%  Set switches for boundary segments to process.

OBC.west  = true;           % process western  boundary segment
OBC.east  = true;           % process eastern  boundary segment
OBC.south = true;           % process southern boundary segment
OBC.north = true;          % process northern boundary segment

S.boundary(1) = OBC.west;
S.boundary(2) = OBC.east;
S.boundary(3) = OBC.south;
S.boundary(4) = OBC.north;

%--------------------------------------------------------------------------
%  Set variables to process.
%--------------------------------------------------------------------------

%  Grid variables.

VarGrd = {'spherical',                                                ...
          'Vtransform', 'Vstretching',                                ...
          'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
          's_rho', 'Cs_r', 's_w', 'Cs_w'};

if (S.spherical),
  if (OBC.west),
    VarGrd = [VarGrd, 'lon_rho_west',  'lat_rho_west',                ...
                      'lon_u_west',    'lat_u_west',                  ...
                      'lon_v_west',    'lat_v_west'];
  end
  if (OBC.east),
    VarGrd = [VarGrd, 'lon_rho_east',  'lat_rho_east',                ...
                      'lon_u_east',    'lat_u_east',                  ...
                      'lon_v_east',    'lat_v_east'];
  end
  if (OBC.south),
    VarGrd = [VarGrd, 'lon_rho_south', 'lat_rho_south',               ...
                      'lon_u_south',   'lat_u_south',                 ...
                      'lon_v_south',   'lat_v_south'];
  end
  if (OBC.north),
    VarGrd = [VarGrd, 'lon_rho_north', 'lat_rho_north',               ...
                      'lon_u_north',   'lat_u_north',                 ...
                      'lon_v_north',   'lat_v_north'];
  end
else
  if (OBC.west),
    VarGrd = [VarGrd, 'x_rho_west',  'y_rho_west',                    ...
                      'x_u_west',    'y_u_west',                      ...
                      'x_v_west',    'y_v_west'];
  end
  if (OBC.east),
    VarGrd = [VarGrd, 'x_rho_east',  'y_rho_east',                    ...
                      'x_u_east',    'y_u_east',                      ...
                      'x_v_east',    'y_v_east'];
  end
  if (OBC.south),
    VarGrd = [VarGrd, 'x_rho_south', 'y_rho_south',                   ...
                      'x_u_south',   'y_u_south',                     ...
                      'x_v_south',   'y_v_south'];
  end
  if (OBC.north),
    VarGrd = [VarGrd, 'x_rho_north', 'y_rho_north',                   ...
                      'x_u_north',   'y_u_north',                     ...
                      'x_v_north',   'y_v_north'];
  end
end

%  ROMS state variables to process.  In 3D applications, the 2D momentum
%  components (ubar,vbar) are compouted by vertically integrating
%  3D momentum component. Therefore, interpolation of (ubar,vbar) is
%  not carried out for efficiency.

VarBry  = {'zeta', 'u', 'v', 'temp', 'salt',...
    'alkalinity','TIC','oxygen','DDCA','CSDC','CLDC','SDOC','SDON','LDOC','LDON',...
    'DDSi','DD_C','DD_N','BAC_','Z2_C','Z2_N','Z1_C','Z1_N','S3CH','S3_C','S3_N',...
    'S2CH','S2_C','S2_N','S1CH','S1_C','S1_N','PO4','SiOH4','NH4','NO3'};
VarList = [VarBry, 'ubar', 'vbar'];

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


%  Get Target grid structure, T.

T = get_roms_grid(GRDname, S);

%%% add by zsm

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
%%%% end of zsm


%  If vector rotation is required in the parent grid, interpolate
%  rotation angle (parent to target) and add it to target grid
%  structure.

T.parent_angle = roms2roms(NWAgrid, P, T, 'angle', [], Rvector,       ...
                           method, offset, RemoveNaN);

%---------------------------------------------------------------------------
%  Create boundary condition Netcdf file.
%---------------------------------------------------------------------------

if (CREATE),

  [status]=c_boundary_cosine(S);

%  Set attributes for "bry_time".

  avalue='seconds since 2012-06-09 01:00:00';
  [status]=nc_attadd(BRYname,'units',avalue,'bry_time');
  
  avalue='gregorian';
  [status]=nc_attadd(BRYname,'calendar',avalue,'bry_time');

%  Set global attribute.

  avalue='Northern Gulf of Mexico, ~2.5 km resolution, Grid a';
  [status]=nc_attadd(BRYname,'title',avalue);

  avalue='ROMS NWA application';
  [status]=nc_attadd(BRYname,'source',avalue);

  [status]=nc_attadd(BRYname,'grd_file',GRDname);
  
%  Write out grid data.
  
  for var = VarGrd,
    field = char(var);
    [err.(field)] = nc_write(BRYname, field, T.(field));
  end

end,
ncdisp(BRYname)
%---------------------------------------------------------------------------
%  Interpolate boundary conditions from Mercator data to application grid.
%---------------------------------------------------------------------------

disp(' ');
disp(['**************************************************************']);
disp(['** Interpolating boundaries conditions from NWA to GOM grid **']);
disp(['**************************************************************']);

%  The NWA data has a time coordinate in seconds, which starts
%  on 1-Jan-1900.

time = nc_read(NWAdata,'ocean_time');
epoch = datenum('01-Jan-1900 00:00');

%  Set reference time such that boundary conditions time are in
%  "seconds since 2000-01-01 00:00:00".

ref_time =datenum('01-Jan-1900');

%  Determine dataset time record to process.

StrRec = find(time == (StrDay-epoch)*86400);
EndRec = find(time == (EndDay-epoch)*86400);

%%%%%%%% 识别不了第一个小时？change by zsm
% if (isempty(StrRec)),
%   error(['  Unable to find starting date in dataset: ', datestr(Tstr)]);
% end
% end
if (isempty(EndRec)),
  error(['  Unable to find ending date in dataset: ', datestr(Tend)]);
end

%%%changeby zsm
% for Rec = StrRec:EndRec,
for Rec = 1:EndRec,

% end
  mydate = datestr(epoch+time(Rec)/86400);

  disp(' ')
  disp(['** Processing: ',mydate,' **']);

%  Interpolate boundary conditions.

  B = obc_roms2roms_zsm(NWAdata, P, T, VarBry, Rec, OBC,                  ...
                    method, offset, RemoveNaN);

%  Extract 3D momentum from boundary structure, B.  

  for var = {'west','east','south','north'},
    edge = char(var);
    ufield = strcat('u','_',edge);
    vfield = strcat('v','_',edge);
    if (OBC.(edge)),
      u.(edge) = B.(ufield);
      v.(edge) = B.(vfield);
    end
  end

%  Set Target Grid vertical level thicknesses, Hz.  We have the
%  choice of using time dependent values (zeta ~= 0) or unperturbed
%  depths (zeta = 0).
  
  if (TimeDependent),
    zeta = roms2roms(NWAdata, P, T, 'zeta', Rec, Rvector,             ...
                     method, offset, RemoveNaN);

    N = S.N;
    igrid = 5;
    [z_w] = set_depth(T.Vtransform, T.Vstretching,                    ...
                      T.theta_s, T.theta_b, T.hc, N,                  ...
                      igrid, T.h, zeta, report);

    Hz = z_w(:,:,2:N+1) -z_w(:,:,1:N);
  else
    Hz = T.Hz;
  end

%  Vertically integrate 3D momentum to compute (ubar,vbar).
  
  [ubar,vbar] = uv_barotropic(u, v, Hz, S.boundary);

%  Load lateral boundary 2D momentum (ubar,vbar) to structure.

  for var = {'west','east','south','north'},
    edge = char(var);
    ufield = strcat('ubar','_',edge);
    vfield = strcat('vbar','_',edge);
    if (OBC.(edge)),
      B.(ufield) = ubar.(edge);
      B.(vfield) = vbar.(edge);
    end
  end

%  Set boundary conditions time (seconds since 2019-01-01 00:00:00).

  % B.bry_time = time(Rec) - ref_time;
  B.bry_time = time(Rec)-(REFTIME).*86400;
  
%  Write out boundary conditions.

  if (WRITE),

    disp(' ');

    BryRec = BryRec+1;

    varlist = fieldnames(B)';
    
    for var = varlist,
      field = char(var);
      [err.(field)] = nc_write(BRYname, field, B.(field), BryRec);
    end

  end

%  Process next boundary record. If processing OpenDAP files, force Java
%  garbage collection.

  [~,url,~] = nc_interface(NWAdata);
  if (url),
    java.lang.System.gc
  end

end