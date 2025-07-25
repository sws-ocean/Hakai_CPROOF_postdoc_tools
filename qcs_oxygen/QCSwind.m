clear

%% 2022

fnameU = 'uwnd.10m.2022.nc';     % 10-m east-wind
fnameV = 'vwnd.10m.2022.nc';     % 10-m north-wind  (same folder)

% 1. Identify the grid point closest to your target location
lat0 = 51;          % °N  (≈ mid-shelf-break in QCS; edit to taste)
lon0 = -130;        % °E  (west is negative)

lat = ncread(fnameU,'lat');      % size [349 277]
lon = ncread(fnameU,'lon');

d2        = (lat - lat0).^2 + (lon - lon0).^2;     % squared distance
[~,idx1d] = min(d2(:));                            % nearest cell
[ix,iy]   = ind2sub(size(lat),idx1d);              % 2-D indices

% 2. Read the full 2022 time series at that single grid cell
u = squeeze(ncread(fnameU,'uwnd' ,[ix iy 1] ,[1 1 Inf]));   % eastward m s⁻¹
v = squeeze(ncread(fnameV,'vwnd' ,[ix iy 1] ,[1 1 Inf]));   % northward m s⁻¹

% 3. Convert model time to MATLAB datetime
time_raw = ncread(fnameU,'time');                 % hours since 1800-01-01
t        = datetime(1800,1,1) + hours(time_raw);  % 3-hourly stamps in 2022

w22.speed = hypot(u,v);                   % m s⁻¹
w22.dir   = mod(atan2d(u,v) + 360,360);   % meteorological ° blowing *from*
w22.lon=lon(ix,iy);
w22.lat=lat(ix,iy);
w22.t=t;

%% 2023
fnameU = 'uwnd.10m.2023.nc';     % 10-m east-wind
fnameV = 'vwnd.10m.2023.nc';     % 10-m north-wind  (same folder)

% 1. Identify the grid point closest to your target location
lat0 = 51.50;          % °N  (≈ mid-shelf-break in QCS; edit to taste)
lon0 = -129.30;        % °E  (west is negative)

lat = ncread(fnameU,'lat');      % size [349 277]
lon = ncread(fnameU,'lon');

d2        = (lat - lat0).^2 + (lon - lon0).^2;     % squared distance
[~,idx1d] = min(d2(:));                            % nearest cell
[ix,iy]   = ind2sub(size(lat),idx1d);              % 2-D indices

% 2. Read the full 2022 time series at that single grid cell
u = squeeze(ncread(fnameU,'uwnd' ,[ix iy 1] ,[1 1 Inf]));   % eastward m s⁻¹
v = squeeze(ncread(fnameV,'vwnd' ,[ix iy 1] ,[1 1 Inf]));   % northward m s⁻¹

% 3. Convert model time to MATLAB datetime
time_raw = ncread(fnameU,'time');                 % hours since 1800-01-01
t        = datetime(1800,1,1) + hours(time_raw);  % 3-hourly stamps in 2022

w23.speed = hypot(u,v);                   % m s⁻¹
w23.dir   = mod(atan2d(u,v) + 360,360);   % meteorological ° blowing *from*
w23.lon=lon(ix,iy);
w23.lat=lat(ix,iy);
w23.t=t;
save('/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/code/R1/QCSwind.mat', ...
    'w22','w23');