%% Think about bifurcation index

% addpath(genpath('/Users/samst/Dropbox/Hakai/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
% addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));
%%
clear

% filename1 = 'UV1cmems_mod_glo_phy_my_0.083deg_P1M-m_1745419519458.nc';
% filename2 = 'UV2cmems_mod_glo_phy_my_0.083deg_P1M-m_1745419691643.nc';
% filename3 = 'UV3cmems_mod_glo_phy_myint_0.083deg_P1M-m_1745419860853.nc';

filename1 = 'Zcmems_mod_glo_phy_my_0.083deg_P1M-m_1746782676267.nc';
filename2 = 'Zcmems_mod_glo_phy_my_0.083deg_P1M-m_1746782774883.nc';
filename3 = 'Zcmems_mod_glo_phy_myint_0.083deg_P1M-m_1746782464799.nc';

% Get info about variables
info = ncinfo(filename1);
vars = {info.Variables.Name};

% Initialize structure
gl = struct();

% Loop through all variables
for i = 1:length(vars)
    varname = vars{i};
    
    % Read data from both files
    data1 = ncread(filename1, varname);
    data2 = ncread(filename2, varname);
    data3 = ncread(filename3, varname);

    % Get dimension names
    dims = info.Variables(i).Dimensions;
    dim_names = string({dims.Name});
    
    % Concatenate along 4th dimension if 'time' is present as 4th dim
    if any(dim_names == "time")
        if length(dims)==1
            gl.(varname) = cat(1, data1, data2, data3);
        elseif length(dims)==3
            gl.(varname) = cat(3, data1, data2, data3);
        elseif length(dims)==4
            gl.(varname) = cat(4, data1, data2, data3);
        end
    else
        gl.(varname) = data1;
    end
end

% Convert time
gl.time = datenum(datetime(1970,1,1) + seconds(gl.time));

[xx,yy,zz,tt] = ndgrid( ...
    gl.longitude, gl.latitude, gl.depth, gl.time);

SA=gsw_SA_from_SP(gl.so,zz,xx,yy);
CT=gsw_CT_from_pt(SA,gl.thetao);
gl.s_t=gsw_sigma0(SA,CT);

%% Interpolate to σ₀ = 26.6
tic

% 1) ensure double
gl.longitude = double(gl.longitude);
gl.latitude  = double(gl.latitude);
gl.depth     = double(gl.depth);
gl.time      = double(gl.time);

% 2) remove duplicate coords (must do this BEFORE any interpn)
[gl.longitude, idx_lon] = unique(gl.longitude, 'stable');
[gl.latitude,  idx_lat] = unique(gl.latitude,  'stable');
[gl.depth,     idx_dep] = unique(gl.depth,     'stable');
[gl.time,      idx_tim] = unique(gl.time,      'stable');

% slice your 4-D fields to match
gl.thetao = gl.thetao(idx_lon, idx_lat, idx_dep, idx_tim);
gl.s_t  = gl.s_t(idx_lon, idx_lat, idx_dep, idx_tim);

% 3) sizes after dedupe
nx = numel(gl.longitude);
ny = numel(gl.latitude);
nz = numel(gl.depth);
nt = numel(gl.time);

% 4) build σ₀ vs. depth matrix [nz × (nx·ny·nt)]
s2 = permute(gl.s_t, [3,1,2,4]);      % → [nz × nx × ny × nt]
s2 = reshape(s2, nz, []);            % → [nz × N],  N = nx·ny·nt

% 5) invert σ₀→depth by looping only where σ₀ brackets 26.6
vq    = 26.6;
N     = size(s2,2);
zv    = nan(1, N);
dmin  = min(s2,[],1);
dmax  = max(s2,[],1);
cols  = find(dmin <= vq & dmax >= vq);

for ii = cols
    sigma_col = s2(:,ii);
    valid     = ~isnan(sigma_col);
    zv(ii)    = interp1( ...
        sigma_col(valid), ...    % X: σ₀ profile
        gl.depth(valid), ...      % V: depth levels
        vq, ...                   % Xq
        'linear', 'extrap');
end

% 6) reshape to 3-D [lon × lat × time]
z_iso = reshape(zv, nx, ny, nt);

% 7) build query grid & bulk‐interpolate u and v
[Li, La, Ti] = ndgrid(gl.longitude, gl.latitude, gl.time);

% uo_iso = interpn( ...
%     gl.longitude, gl.latitude, gl.depth, gl.time, gl.uo, ...
%     Li, La, z_iso, Ti, ...
%     'linear', NaN);
% 
% vo_iso = interpn( ...
%     gl.longitude, gl.latitude, gl.depth, gl.time, gl.vo, ...
%     Li, La, z_iso, Ti, ...
%     'linear', NaN);

% T_iso = interpn( ...
%     gl.longitude, gl.latitude, gl.depth, gl.time, gl.thetao, ...
%     Li, La, z_iso, Ti, ...
%     'linear', 'extrap');

% build the interpolant with linear extrapolation
Ftheta = griddedInterpolant( ...
    {gl.longitude, gl.latitude, gl.depth, gl.time}, ...
    gl.thetao, ...
    'linear', ...    % interpolation method
    'linear' );      % extrapolation method

% evaluate on the isopycnal grid
T_iso = Ftheta( Li, La, z_iso, Ti );

% 8) diagnostics
nTotal = numel(z_iso);
nIso   = nnz(~isnan(z_iso));
nNaN   = nnz( isnan(z_iso));

fprintf('Total grid pts: %d\n',        nTotal);
fprintf('Iso crossings : %d (%.1f%%)\n', nIso,  100*nIso/nTotal);
fprintf('NaN points    : %d (%.1f%%)\n', nNaN, 100*nNaN/nTotal);

toc


%% Make some plots to stare at
yr=year(gl.time);
yrs=yr(1):2024;

col=[255 214 140]/255; % YELLOW
lon_lim=[min(gl.longitude(:)) max(gl.longitude(:))];
lat_lim=[min(gl.latitude) max(gl.latitude)];
[xx,yy]=meshgrid(gl.longitude,gl.latitude);

gl.aTo=NaN([length(gl.latitude) length(gl.longitude) length(yrs)]);

aMtime=NaN(1,length(yrs));
for i=1:length(yrs)
    msk=yr==yrs(i);
    tmpT=mean(squeeze(T_iso(:,:,msk)),3,'omitnan');
    gl.aTo(:,:,i)=tmpT'-nanmean(tmpT(:));
    aMtime(i)=datenum(yrs(i),06,01);

    % f1=figure('units','centimeters','outerposition',...
    %     [0.01 0.01 8 8],'color','w');
    % m_proj('mercator', 'lon', lon_lim, 'lat', lat_lim); % North Pacific region
    % % text(0.025, 0.95, label{i}, 'Units', 'normalized');
    % m_contourf(gl.longitude,gl.latitude,tmpT'-nanmean(tmpT(:)),-4:0.1:4,...
    %     'linestyle','none');
    % clim([-1 1]);
    % colormap(cmocean('balance','pivot',0))
    % m_gshhs_i('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
    % title(sprintf('%4.0f: annual mean temperature anomaly (T'') at 400m',yrs(i)));
    % c=colorbar;
    % c.Label.String = 'T'' (m s^{-1})';
    % m_grid('tickdir', 'out','linestyle','none'); % Add grid and labels
    % set(findall(gcf,'-property','fontsize'),'fontsize',8);
    % set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
end



%% DO PCA
%— assume gl.aTo(i,j,t), gl.longitude (1×nj), gl.latitude (1×ni) in workspace
[ni,nj,nt] = size(gl.aTo);

% 1) flatten and mask NaNs
Xall = reshape(gl.aTo, ni*nj, nt);        
good = all(~isnan(Xall), 2);                
X = Xall(good, :);

% 2) remove temporal mean
mT = mean(X, 2);
X = X - mT;

% 3) prep for PCA: rows=time, cols=space
Xp = X.';   % nt × nGood

% 4) EOF/PCA
[coeff, score, ~, ~, explained] = pca(Xp);

% 5) reconstruct EOFs on full grid
nModes = min(4, size(coeff,2));
EOFall = nan(ni*nj, nModes);
EOFall(good, :) = coeff(:, 1:nModes);
EOF = reshape(EOFall, ni, nj, nModes);

% 6) prepare lat/lon grids
[Lon, Lat] = meshgrid(gl.longitude, gl.latitude);

% 7) set limits & coast color
lon_lim = [min(gl.longitude) max(gl.longitude)];
lat_lim = [min(gl.latitude)  max(gl.latitude)];
coast_col = [0.7 0.7 0.7];

% 8) plot EOF patterns with m_map
for k = 1:nModes
    figure
    m_proj('mercator', 'lon', lon_lim, 'lat', lat_lim);
    m_grid('tickdir','out','linestyle','none');

    % anomaly about its own mean
    Z = EOF(:,:,k);

    % filled contour on map
    m_contourf(Lon, Lat, Z, 'linestyle','none');

    colormap( m_colmap('jet'));
    m_coast('patch', coast_col);

    title(sprintf('EOF %d (%.1f%% var)', k, explained(k)));
    colorbar('southoutside')
end

% 9) plot PC time series
figure
plot(1:nt, score(:,1:nModes), 'LineWidth', 1.2)
legend(arrayfun(@(x) sprintf('PC%d', x), 1:nModes, 'uni',0))
xlabel('Time index'), ylabel('Amplitude')
title('Principal Component time series')

% save('/Users/samuelstevens/Dropbox/Hakai/gyre/supply/data/NEPGpcsISO.mat','score');

%% What's a physical representation of PC2? 

%— after you have EOF (ni×nj×nModes) and score (nt×nModes) from your PCA

% 1) compute the climatological mean temperature field
meanField = squeeze(nanmean(gl.aTo, 3));  % ni×nj

% 2) extract EOF2 spatial pattern
EOF2 = EOF(:,:,1);                          % ni×nj

% 3) compute PC2 series and its standard deviation
PC2 = score(:,1);                          
sd2 = std(PC2);

% 4) define amplitudes for “strongly positive”, “neutral”, “strongly negative”
amp   = 2*sd2;                              % choose 2σ for “strong”
mults = [ amp, -amp];                    
labels = {'+2σ','-2σ'};

% 5) build 2-D lon/lat grids (if not already)
[Lon, Lat] = meshgrid(gl.longitude, gl.latitude);

% 6) plot the three maps
lon_lim = [min(gl.longitude) max(gl.longitude)];
lat_lim = [min(gl.latitude)  max(gl.latitude)];

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 19 21],'color','w');

label={'a)';'b)';'c)';'d)';'e)';'f)'};
t = tiledlayout(2,3,'TileSpacing', 'Compact', 'Padding', 'Compact');
nexttile;
m_proj('mercator', 'lon', lon_lim, 'lat', lat_lim);

% anomaly about its own mean
Z = EOF(:,:,1);

% filled contour on map
m_contourf(Lon, Lat, Z,'levelstep',0.0005, 'linestyle','none');

colormap(gca,cmocean('curl','pivot',0));
m_gshhs_h('patch',[0.5 0.5 0.5],'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('tickdir','out','linestyle','none');

title(sprintf('EOF %d (%.1f%% var)', 2, explained(2)));
c=colorbar('southoutside');
c.Label.String = 'EOF2 amplitude (°C)';
text(0.925, 0.95, label{1}, 'Units', 'normalized');

for m = 1:length(mults)
    nexttile;
    % reconstruct the temperature field at this PC amplitude
    Tmap = meanField + mults(m).*EOF2;

    m_proj('mercator','lon',lon_lim,'lat',lat_lim);

    % choose contour levels around your data range
    cli = [min(Tmap(:)), max(Tmap(:))];
    levels = linspace(cli(1), cli(2), 40);

    m_contourf(Lon, Lat, Tmap, levels, 'linestyle','none');
    m_gshhs_h('patch',[0.5 0.5 0.5],'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
    m_grid('tickdir','out','linestyle','none','yticklabels',[]);
    title(sprintf('PC2 = %s', labels{m}));
    caxis([-1 1])
    c=colorbar('southoutside');
    c.Label.String = 'T'' (°C)';
    colormap(gca,cmocean('balance','pivot',0));
    text(0.925, 0.95, label{m+1}, 'Units', 'normalized');
end

ax1=axes('position',[0.2 0.3 0.6 0.15]); hold on
plot(yrs, score(:,2),'k','LineWidth', 1.2);
scatter(yrs, score(:,2),10,'filled','k');
plot(yrs, score(:,[1 3 4]),'color',[.8 .8 .8 .5],'LineWidth', 1);
xlabel('Time'), ylabel('PC2');
grid on; axis tight
text(0.975, 0.95, label{4}, 'Units', 'normalized');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/EOFV2.pdf -dpdf -nofontswap

% %% compare pcs with BI
% 
% BI=load('BifurcationIndexApr25.mat');
% for i=1:length(BI.BI)
%     BI.mtime(i)=datenum(BI.yr(i),1,1);
% end
% 
% figure;
% subplot(3,1,1)
% plot(BI.yr,BI.BI)
% yyaxis right
% plot(yrs,score(:,1))
% grid on;
% subplot(3,1,2)
% plot(BI.yr,BI.BI)
% yyaxis right
% plot(yrs,score(:,2))
% grid on;
% subplot(3,1,3)
% plot(BI.yr,BI.BI)
% yyaxis right
% plot(yrs,score(:,3))
% grid on;
% 

%% Look at monthly PCA
%— assume gl.aTo(i,j,t), gl.longitude (1×nj), gl.latitude (1×ni) in workspace
% mnthTo=permute(squeeze(gl.thetao),[2,1,3]);
mnthTo=permute(squeeze(T_iso),[2,1,3]);
[ni,nj,nt] = size(mnthTo);

% 1) flatten and mask NaNs
Xall = reshape(mnthTo, ni*nj, nt);    

% Interp across months

for i=1:size(Xall)
    Xall(i,:)=inpaint_nans(Xall(i,:),5);
end

good = all(~isnan(Xall), 2);                
X = Xall(good, :);

% 2) remove temporal mean
mT = mean(X, 2);
X = X - mT;

% 3) prep for PCA: rows=time, cols=space
Xp = X.';   % nt × nGood

% 4) EOF/PCA
[coeff, score, ~, ~, explained] = pca(Xp);

% 5) reconstruct EOFs on full grid
nModes = min(4, size(coeff,2));
EOFall = nan(ni*nj, nModes);
EOFall(good, :) = coeff(:, 1:nModes);
EOF = reshape(EOFall, ni, nj, nModes);

% 6) prepare lat/lon grids
[Lon, Lat] = meshgrid(gl.longitude, gl.latitude);

% 7) set limits & coast color
lon_lim = [min(gl.longitude) max(gl.longitude)];
lat_lim = [min(gl.latitude)  max(gl.latitude)];
coast_col = [0.7 0.7 0.7];

% 8) plot EOF patterns with m_map
for k = 1:nModes
    figure
    m_proj('mercator', 'lon', lon_lim, 'lat', lat_lim);
    m_grid('tickdir','out','linestyle','none');

    % anomaly about its own mean
    Z = EOF(:,:,k);

    % filled contour on map
    m_contourf(Lon, Lat, Z,'levelstep',0.0005, 'linestyle','none');

    colormap( m_colmap('jet'));
    m_coast('patch', coast_col);

    title(sprintf('EOF %d (%.1f%% var)', k, explained(k)));
    colorbar('southoutside')
end

% 9) plot PC time series
figure
plot(1:nt, score(:,1:nModes), 'LineWidth', 1.2)
legend(arrayfun(@(x) sprintf('PC%d', x), 1:nModes, 'uni',0))
xlabel('Time index'), ylabel('Amplitude')
title('Principal Component time series')

aScore=NaN(32,386);
yrs=year(gl.time);
for i=1993:2024
    msk=yrs==i;
    aScore(i-1992,:)=mean(score(msk,:),1);
end

% save('/Users/samuelstevens/Dropbox/Hakai/gyre/supply/data/NEPGpcsMONTH.mat','score');


%% What's a physical representation of PC2? 

%— after you have EOF (ni×nj×nModes) and score (nt×nModes) from your PCA

% 1) compute the climatological mean temperature field
meanField = squeeze(nanmean(gl.aTo, 3));  % ni×nj

EOFidx=3;

% 2) extract EOF2 spatial pattern
EOF2 = EOF(:,:,EOFidx);                          % ni×nj

% 3) compute PC2 series and its standard deviation
PC2 = score(:,EOFidx);                          
sd2 = std(PC2);

% 4) define amplitudes for “strongly positive”, “neutral”, “strongly negative”
amp   = 2*sd2;                              % choose 2σ for “strong”
mults = [ amp, -amp];                    
labels = {'+2σ','-2σ'};

% 5) build 2-D lon/lat grids (if not already)
[Lon, Lat] = meshgrid(gl.longitude, gl.latitude);

% 6) plot the three maps
lon_lim = [min(gl.longitude) max(gl.longitude)];
lat_lim = [min(gl.latitude)  max(gl.latitude)];

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 19 21],'color','w');

label={'a)';'b)';'c)';'d)';'e)';'f)'};
t = tiledlayout(2,3,'TileSpacing', 'Compact', 'Padding', 'Compact');
nexttile;
m_proj('mercator', 'lon', lon_lim, 'lat', lat_lim);

% anomaly about its own mean
Z = EOF(:,:,EOFidx);

% filled contour on map
m_contourf(Lon, Lat, Z,'levelstep',0.0005, 'linestyle','none');

colormap(gca,cmocean('curl','pivot',0));
m_gshhs_h('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('tickdir','out','linestyle','none');

title(sprintf('EOF %d (%.1f%% var)', EOFidx, explained(EOFidx)));
c=colorbar('southoutside');
c.Label.String = 'EOF2 amplitude (°C)';
text(0.925, 0.95, label{1}, 'Units', 'normalized');

for m = 1:length(mults)
    nexttile;
    % reconstruct the temperature field at this PC amplitude
    Tmap = meanField + mults(m).*EOF2;

    m_proj('mercator','lon',lon_lim,'lat',lat_lim);

    % choose contour levels around your data range
    cli = [min(Tmap(:)), max(Tmap(:))];
    levels = linspace(cli(1), cli(2), 40);

    m_contourf(Lon, Lat, Tmap, levels, 'linestyle','none');
    m_gshhs_h('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
    m_grid('tickdir','out','linestyle','none','yticklabels',[]);
    title(sprintf('PC%.0f %s',EOFidx, labels{m}));
    caxis([-1 1])
    c=colorbar('southoutside');
    c.Label.String = 'T'' (°C)';
    colormap(gca,cmocean('balance','pivot',0));
    text(0.925, 0.95, label{m+1}, 'Units', 'normalized');
end

scoreM=[];
for i=1:length(yrs)
    msk=yr==yrs(i);
    scoreM(i,:)=nanmean(score(msk,:));
end

ax1=axes('position',[0.2 0.3 0.6 0.15]); hold on
plot(gl.time, score(:,EOFidx),'color',[0 0 0 0.5],'LineWidth', 0.8);
% plot(aMtime, scoreM(:,EOFidx),'k','LineWidth', 1.2);
plot(gl.time, smoothdata(score(:,EOFidx),'movmean',12),'k','LineWidth', 1.2)
xlabel('Time'), ylabel('PC4');
grid on; axis tight; axdate(20);
text(0.975, 0.95, label{4}, 'Units', 'normalized');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/EOF_Monthly.pdf -dpdf -nofontswap

%% Compare to other things

% figure; hold on;
% plot(gl.time,score(:,3));
% plot(aMtime,smoothdata(scoreM(:,3),'movmean',1),'linewidth',2);
% a=load('NEPGpcsM.mat');
% plot(a.aMtime,a.scoreM(:,3)','linewidth',2,'color','g');
% 
% figure;
% scatter(score(:,4),a.score(:,3),'k');
% hold on
% scatter(scoreM(:,4),a.scoreM(:,3),'r','filled');
% 
% clc
% for i=1:5
%     [r1,p1]=corrcoef(score(:,i),a.score(:,3));
%     fprintf('Mode %.0f: r=%.2f, p=%.2f\n',i,r1(2),p1(2));
% end

clc
load allSpice.mat
for i=1:5
    msk=~isnan(P4.aSpice(2,:));
    [r1,p1]=corrcoef(aScore(msk,i),P4.aSpice(2,msk));
    fprintf('Mode %.0f vs P4 spice: r=%.2f, p=%.2f\n',i,r1(2),p1(2));
end
disp(' ');

for i=1:5
    msk=~isnan(CS01.aSpice(2,:));
    [r1,p1]=corrcoef(aScore(msk,i),CS01.aSpice(2,msk));
    fprintf('Mode %.0f vs CS01 spice: r=%.2f, p=%.2f\n',i,r1(2),p1(2));
end
disp(' ');

for i=1:5
    msk=~isnan(CS09.aSpice(2,:));
    [r1,p1]=corrcoef(aScore(msk,i),CS09.aSpice(2,msk));
    fprintf('Mode %.0f vs CS09 spice: r=%.2f, p=%.2f\n',i,r1(2),p1(2));
end


disp(' ');
load breathBifurcate.mat
for i=1:5
    [r1,p1]=corrcoef(score(:,i),mode_time_series(:,1));
    fprintf('Mode %.0f vs Breathing Mode spice: r=%.2f, p=%.2f\n',i,r1(2),p1(2));
end

disp(' ');
for i=1:5
    [r1,p1]=corrcoef(score(:,i),mode_time_series(:,2));
    fprintf('Mode %.0f vs Bifurcation Mode spice: r=%.2f, p=%.2f\n',i,r1(2),p1(2));
end

%%
save('/Users/samuelstevens/Dropbox/Hakai/gyre/supply/data/NEPGpcsMISO.mat','score','scoreM','aMtime');


