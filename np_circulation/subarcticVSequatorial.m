%% Recreate TnK

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

%%
clear

%% Load in PEW data

% Initialize structure
G= struct();

info = ncinfo('PEW1.nc');
vars = {info.Variables.Name};

s={'PEW'};

% Loop through all variables
for f=1:length(s)
    for i=1:length(vars)
        varname = vars{i};

        % Read data from both files
        data1 = ncread([s{f},'1.nc'],varname);
        data2= ncread([s{f},'2.nc'],varname);
        data3= ncread([s{f},'3.nc'],varname);
        data4= ncread([s{f},'4.nc'],varname);

        % Get dimension names
        dims = info.Variables(i).Dimensions;
        dim_names = string({dims.Name});

        % Concatenate along 4th dimension if 'time' is present as 4th dim
        if any(dim_names == "time")
            if length(dims)==1
                eval([s{f},'.(varname) = double(cat(1, data1, data2, data3, data4));'])
            elseif length(dims)==3
                eval([s{f},'.(varname) = double(squeeze(cat(3, data1, data2, data3, data4)));'])
            elseif length(dims)==4
                eval([s{f},'.(varname) = double(squeeze(cat(4, data1, data2, data3, data4)));'])
            end
        else
                eval([s{f},'.(varname) = double(data1);'])
        end
    end
end

PEWt=squeeze(mean(PEW.thetao,[1 2 4],'omitnan'));
PEWs=squeeze(mean(PEW.so,[1 2 4],'omitnan'));
PEWz=PEW.depth;

clearvars -except PEWs PEWt PEWz

%% Load in SAW data

% Initialize structure
G= struct();

info = ncinfo('SAW1.nc');
vars = {info.Variables.Name};

s={'SAW'};

% Loop through all variables
for f=1:length(s)
    for i=1:length(vars)
        varname = vars{i};

        % Read data from both files
        data1 = ncread([s{f},'1.nc'],varname);
        data2= ncread([s{f},'2.nc'],varname);
        data3= ncread([s{f},'3.nc'],varname);
        data4= ncread([s{f},'4.nc'],varname);

        % Get dimension names
        dims = info.Variables(i).Dimensions;
        dim_names = string({dims.Name});

        % Concatenate along 4th dimension if 'time' is present as 4th dim
        if any(dim_names == "time")
            if length(dims)==1
                eval([s{f},'.(varname) = double(cat(1, data1, data2, data3, data4));'])
            elseif length(dims)==3
                eval([s{f},'.(varname) = double(squeeze(cat(3, data1, data2, data3, data4)));'])
            elseif length(dims)==4
                eval([s{f},'.(varname) = double(squeeze(cat(4, data1, data2, data3, data4)));'])
            end
        else
                eval([s{f},'.(varname) = double(data1);'])
        end
    end
end

SAWt=squeeze(mean(SAW.thetao,[1 2 4],'omitnan'));
SAWs=squeeze(mean(SAW.so,[1 2 4],'omitnan'));
SAWz=SAW.depth;

clearvars -except SAWs SAWt SAWz PEWs PEWt PEWz

%% Interpolate onto regular grids

SA=gsw_SA_from_SP(SAWs,SAWz,ones(size(SAWt)).*-142.5,ones(size(SAWt)).*50);
CT=gsw_CT_from_pt(SA,SAWt);
SAW_st=gsw_sigma0(SA,CT);

SA=gsw_SA_from_SP(PEWs,PEWz,ones(size(PEWt)).*-142.5,ones(size(PEWt)).*50);
CT=gsw_CT_from_pt(SA,PEWt);
PEW_st=gsw_sigma0(SA,CT);

sgrid=[23:0.05:27.5]';

msk=~isnan(SAW_st);
SAWt=interp1(SAW_st(msk),SAWt(msk),sgrid);
SAWs=interp1(SAW_st(msk),SAWs(msk),sgrid);
SAWz=interp1(SAW_st(msk),SAWz(msk),sgrid);

msk=~isnan(PEW_st);
PEWt=interp1(PEW_st(msk),PEWt(msk),sgrid);
PEWs=interp1(PEW_st(msk),PEWs(msk),sgrid);
PEWz=interp1(PEW_st(msk),PEWz(msk),sgrid);


%% Create mixtures
mixGridT=NaN(length(sgrid),9);
mixGridS=NaN(length(sgrid),9);

for i=1:9
    mixGridT(:,i)=SAWt.*(i/10)+PEWt.*(1-(i/10));
    mixGridS(:,i)=SAWs.*(i/10)+PEWs.*(1-(i/10));
end


%% figure
load allSpice.mat

% Calculate spiciness and density grid
Sgrid = 32:0.01:35.5;
Tgrid = 2:0.1:20;
[SS, TT] = meshgrid(Sgrid, Tgrid);

SA_grid = gsw_SA_from_SP(SS, 0, -142.5, 50);  % surface pressure
CT_grid = gsw_CT_from_t(SA_grid, TT, 0);
sigma0 = gsw_sigma0(SA_grid, CT_grid);
spice = gsw_spiciness0(SA_grid, CT_grid);

% Plot setup
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 15 21],'color','w');
hold on;

% Background spice and sigma contours
[h1, c1] = contour(SS, TT, spice, 'k:');  % spiciness
[h2,c2] = contour(SS, TT, sigma0, 'k-'); % sigma0

% Plot endmember profiles
l1=plot(SAWs, SAWt, 'b', 'LineWidth', 2);
l2=plot(PEWs, PEWt, 'r', 'LineWidth', 2);

% plot CS01 and P4
scatter(CS01.sal,CS01.temp,1,'filled','markerfacecolor', ...
    rgb_x('rose'),'markerfacealpha',0.2);
scatter(P4.sal,P4.temp,1,'filled','markerfacecolor', ...
    rgb_x('steel blue'),'markerfacealpha',0.2);
plot(mean(P4.sal,2,'omitnan'),mean(P4.temp,2,'omitnan'),'w','linewidth',5);
l3=plot(mean(P4.sal,2,'omitnan'),mean(P4.temp,2,'omitnan'),'linewidth',3, ...
    'color',rgb_x('steel blue'));
plot(mean(CS01.sal,2,'omitnan'),mean(CS01.temp,2,'omitnan'),'w','linewidth',5);
l4=plot(mean(CS01.sal,2,'omitnan'),mean(CS01.temp,2,'omitnan'),'linewidth',3, ...
    'color',rgb_x('rose'));

% Plot mixing lines
for i = 1:size(mixGridT,2)
    plot(mixGridS(:,i), mixGridT(:,i),'--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
end

% Style
xlabel('Salinity (psu)');
ylabel('Temperature (°C)');
% title('T–S Diagram with Spiciness and σ₀: PEW vs SAW Mixing');
% grid on;
box on;
set(gca, 'FontSize', 10);
xlim([32.5 35.5]);
ylim([2 18]);
 
set(findall(gcf,'-property','fontsize'),'fontsize',10);

clabel(h1, 'FontSize', 8, 'Color', [0.4 0.4 0.4],'fontname','Latin Modern Roman');
clabel(h2, 'FontSize', 8, 'Color', [0 0 0],'fontname','Latin Modern Roman');

legend([c1,c2,l1,l2,l3,l4],{'Spiciness', 'σ₀', 'PSUW', 'PEW','P4','CS01'}, 'Location', 'southwest','FontSize', 8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%

% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/TSmixing.jpg -m4


%% Find fractions using linear mixing model along isopycnals

% CS01
[nz_cs, nprof_cs] = size(CS01.temp);
alpha_SAW_CS01 = nan(nz_cs, nprof_cs);
alpha_PEW_CS01 = nan(nz_cs, nprof_cs);

for i = 1:nz_cs
  for j = 1:nprof_cs
    st_val = CS01.s_t(i,j);
    if isnan(st_val), continue, end

    % interpolate endmember at this isopycnal
    SAW_t = interp1(sgrid, SAWt, st_val);
    SAW_s = interp1(sgrid, SAWs, st_val);
    PEW_t = interp1(sgrid, PEWt, st_val);
    PEW_s = interp1(sgrid, PEWs, st_val);

    % sample T/S
    T_val = CS01.temp(i,j);
    S_val = CS01.sal(i,j);

    % mixing‐vector & sample‐offset
    dT   = SAW_t - PEW_t;
    dS   = SAW_s - PEW_s;
    dT_s = T_val - PEW_t;
    dS_s = S_val - PEW_s;

    % projection for SAW fraction
    alpha = (dT_s*dT + dS_s*dS) / (dT^2 + dS^2);
    alpha = max(min(alpha,1),0);

    alpha_SAW_CS01(i,j) = alpha;
    alpha_PEW_CS01(i,j) = 1 - alpha;
  end
end

% P4
[nz_p4, nprof_p4] = size(P4.temp);
alpha_SAW_P4 = nan(nz_p4, nprof_p4);
alpha_PEW_P4 = nan(nz_p4, nprof_p4);

for i = 1:nz_p4
  for j = 1:nprof_p4
    st_val = P4.s_t(i,j);
    if isnan(st_val), continue, end

    % interpolate endmember at this isopycnal
    SAW_t = interp1(sgrid, SAWt, st_val);
    SAW_s = interp1(sgrid, SAWs, st_val);
    PEW_t = interp1(sgrid, PEWt, st_val);
    PEW_s = interp1(sgrid, PEWs, st_val);

    % sample T/S
    T_val = P4.temp(i,j);
    S_val = P4.sal(i,j);

    % mixing‐vector & sample‐offset
    dT   = SAW_t - PEW_t;
    dS   = SAW_s - PEW_s;
    dT_s = T_val - PEW_t;
    dS_s = S_val - PEW_s;

    % projection for SAW fraction
    alpha = (dT_s*dT + dS_s*dS) / (dT^2 + dS^2);
    alpha = max(min(alpha,1),0);

    alpha_SAW_P4(i,j) = alpha;
    alpha_PEW_P4(i,j) = 1 - alpha;
  end
end


%%  Depth profiles;
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 10 12],'color','w');

% subplot(2,1,1);
plot(alpha_PEW_P4,1:1000,'color',[rgb_x('steel blue') 0.02],'linewidth',0.75);
hold on;
plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'w','linewidth',4);
l1=plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'color',rgb_x('steel blue'), ...
    'linewidth',2);

mnths=month(P4.mtime);
msk=mnths>=4 & mnths<=10;
plot(alpha_PEW_CS01,1:1000,'color',[rgb_x('rose') 0.02],'linewidth',0.75);
plot(mean(alpha_PEW_CS01,2,'omitnan'),1:1000,'w','linewidth',4);
l2=plot(mean(alpha_PEW_CS01,2,'omitnan'),1:1000,'color',rgb_x('rose'), ...
    'linewidth',2);
axis ij;
xlabel('Fraction of PEW');
ylabel('Depth (m)');
grid on;

set(findall(gcf,'-property','fontsize'),'fontsize',10);

% plot density markers
csst=mean(CS01.s_t,2,'omitnan');
p4st=mean(P4.s_t,2,'omitnan');
for i=25:0.2:27
    y=interp1(csst,1:1000,i);
    x=interp1(1:1000,mean(alpha_PEW_CS01,2,'omitnan'),y);
    scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('rose'));

    if i==26.4
        text(x,y,'26.4 \sigma_\theta','fontsize',8);
    end

    y=interp1(p4st,1:1000,i);
    x=interp1(1:1000,mean(alpha_PEW_P4,2,'omitnan'),y);
    s1=scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('steel blue'));

    if i==26.4
        text(x,y,'26.4 \sigma_\theta','fontsize',8);
    end

end

ylim([0 500]);
legend([l1,l2,s1],{'P4','CS01','0.2 \sigma_\theta spacing'}, 'Location', 'southwest','FontSize', 8);


set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


% plot(mean(alpha_SAW_P4(:,msk),2,'omitnan'),1:1000,'w','linewidth',2);
% l=plot(mean(alpha_SAW_P4(:,msk),2,'omitnan'),1:1000,'-','color',rgb_x('steel blue'), ...
%     'linewidth',1);
% % plot(mean(alpha_SAW_P4(:,~msk),2,'omitnan'),1:1000,'w','linewidth',2);
% l=plot(mean(alpha_SAW_P4(:,~msk),2,'omitnan'),1:1000,'-','color',rgb_x('steel blue'), ...
%     'linewidth',1);

% subplot(1,2,2);
% plot(alpha_PEW_P4,1:1000,'color',[rgb_x('steel blue') 0.2],'linewidth',0.75);
% hold on;
% plot(alpha_PEW_CS01,1:1000,'color',[rgb_x('rose') 0.2],'linewidth',0.75);
% 
% plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'w','linewidth',4);
% l=plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'color',rgb_x('steel blue'), ...
%     'linewidth',2);
% plot(mean(alpha_PEW_CS01,2,'omitnan'),1:1000,'w','linewidth',4);
% ll=plot(mean(alpha_PEW_CS01,2,'omitnan'),1:1000,'color',rgb_x('rose'), ...
%     'linewidth',2);
% axis ij;
% xlabel('Fraction of PEW');
% ylabel('Depth (m)');
% grid on;
% legend([l,ll],'P4','CS01','location','best');




%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/SAW_profs.jpg -m4

%% Calculate time variability along isopycnal

% Target isopycnal
ptarget = 26.5;
% CS01: interpolate alpha to σ₀ = 26.5
alpha_SAW_CS01_26p5 = nan(1, nprof_cs);
alpha_PEW_CS01_26p5 = nan(1, nprof_cs);

for j = 1:nprof_cs
  st_prof  = CS01.s_t(:,j);
  saw_prof = alpha_SAW_CS01(:,j);
  pew_prof = alpha_PEW_CS01(:,j);

  valid = ~isnan(st_prof) & ~isnan(saw_prof);
  if sum(valid) >= 2
    % linear interpolation at ptarget, out-of-range → NaN
    alpha_SAW_CS01_26p5(j) = interp1(st_prof(valid), saw_prof(valid), ptarget, 'linear', NaN);
    alpha_PEW_CS01_26p5(j) = interp1(st_prof(valid), pew_prof(valid), ptarget, 'linear', NaN);
  end
end

% P4: interpolate alpha to σ₀ = 26.5
alpha_SAW_P4_26p5 = nan(1, nprof_p4);
alpha_PEW_P4_26p5 = nan(1, nprof_p4);

for j = 1:nprof_p4
  st_prof  = P4.s_t(:,j);
  saw_prof = alpha_SAW_P4(:,j);
  pew_prof = alpha_PEW_P4(:,j);

  valid = ~isnan(st_prof) & ~isnan(saw_prof);
  if sum(valid) >= 2
    alpha_SAW_P4_26p5(j) = interp1(st_prof(valid), saw_prof(valid), ptarget, 'linear', NaN);
    alpha_PEW_P4_26p5(j) = interp1(st_prof(valid), pew_prof(valid), ptarget, 'linear', NaN);
  end
end

%% Plot time variability

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 10 8],'color','w');
hold on;
scatter(CS01.mtime,alpha_SAW_CS01_26p5,20, '^', 'filled','markerfacecolor',...
    rgb_x('rose'),'MarkerFaceAlpha',0.4);
scatter(P4.mtime,alpha_SAW_P4_26p5,20, '^', 'filled','markerfacecolor',...
    rgb_x('steel blue'),'MarkerFaceAlpha',0.2);

hold on;
scatter(CS01.mtime,alpha_PEW_CS01_26p5,20, 'filled','markerfacecolor',...
    rgb_x('rose'),'MarkerFaceAlpha',0.4);
scatter(P4.mtime,alpha_PEW_P4_26p5,20,'filled','markerfacecolor',...
    rgb_x('steel blue'),'MarkerFaceAlpha',0.2);

% Calculate annual means for CS01
[Y_cs,~,~] = datevec(CS01.mtime);
years_cs  = unique(Y_cs);
mean_SAW_cs = nan(size(years_cs));
mean_PEW_cs = nan(size(years_cs));
for k = 1:length(years_cs)
    idx = (Y_cs == years_cs(k));
    mean_SAW_cs(k) = nanmean(alpha_SAW_CS01_26p5(idx));
    mean_PEW_cs(k) = nanmean(alpha_PEW_CS01_26p5(idx));
end

% Plot CS01 annual means
l2=plot(datenum(years_cs,1,1), mean_SAW_cs,'color',rgb_x('rose'),'LineWidth', 1.5);
plot(datenum(years_cs,1,1), mean_PEW_cs,'color',rgb_x('rose'),'LineWidth', 1.5);

% Calculate annual means for P4
[Y_p4,~,~] = datevec(P4.mtime);
years_p4  = unique(Y_p4);
mean_SAW_p4 = nan(size(years_p4));
mean_PEW_p4 = nan(size(years_p4));
for k = 1:length(years_p4)
    idx = (Y_p4 == years_p4(k));
    mean_SAW_p4(k) = nanmean(alpha_SAW_P4_26p5(idx));
    mean_PEW_p4(k) = nanmean(alpha_PEW_P4_26p5(idx));
end

% Plot CS01 annual means
l1=plot(datenum(years_p4,1,1), mean_SAW_p4,'color',rgb_x('steel blue'),'LineWidth', 1.5);
plot(datenum(years_p4,1,1), mean_PEW_p4,'color',rgb_x('steel blue'),'LineWidth', 1.5);

% Labels & formatting
xlabel('Year');
ylabel('Fraction at \sigma_0 = 26.5');
axdate(20);      % format tick labels with 20° rotation
% legend('Location','best');
box on; grid on; axis tight;

set(findall(gcf,'-property','fontsize'),'fontsize',8);
text(0.1,0.5,'PSUW','units','normalized','fontsize',10);
text(0.6,0.5,'PEW','units','normalized','fontsize',10);


legend([l1,l2],{'P4','CS01'}, 'Location', 'southwest','FontSize', 8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/interannual_PEW_SAUW.jpg -m4


%% Does interannual variability correlate with the breathing mode? No
BB = load('breathBifurcate.mat');

yrsP4 = year(P4.mtime);
yrsCS01 = year(CS01.mtime);
yrsBB = year(BB.G.time);
clc

% Initialize
nyrs = 2024 - 1993 + 1;
Ba = NaN(nyrs,1);
CS01a = NaN(nyrs,1);
P4a = NaN(nyrs,1);

for i = 1993:2024
    idx = i - 1992;

    msk = yrsBB == i;
    Ba(idx) = mean(BB.mode_time_series(msk,2), 'omitnan');

    msk = yrsCS01 == i;
    CS01a(idx) = mean(alpha_SAW_CS01_26p5(msk), 'omitnan');

    msk = yrsP4 == i;
    P4a(idx) = mean(alpha_SAW_P4_26p5(msk), 'omitnan');
end

% Standard (contemporaneous) correlation
disp('Contemporaneous correlations:')
msk = ~isnan(CS01a) & ~isnan(Ba);
[r,p] = corrcoef(CS01a(msk), Ba(msk));
disp(['CS01: r = ', num2str(r(1,2)), ', p = ', num2str(p(1,2))])

msk = ~isnan(P4a) & ~isnan(Ba);
[r,p] = corrcoef(P4a(msk), Ba(msk));
disp(['P4:   r = ', num2str(r(1,2)), ', p = ', num2str(p(1,2))])

% Lagged correlation (Ba leads by 1 year)
disp('Lagged correlations (Ba leads by 1 year):')
Ba_lag1 = Ba(1:end-1);
CS01a_lag = CS01a(2:end);
P4a_lag = P4a(2:end);

msk = ~isnan(Ba_lag1) & ~isnan(CS01a_lag);
[r,p] = corrcoef(Ba_lag1(msk), CS01a_lag(msk));
disp(['CS01: r = ', num2str(r(1,2)), ', p = ', num2str(p(1,2))])

msk = ~isnan(Ba_lag1) & ~isnan(P4a_lag);
[r,p] = corrcoef(Ba_lag1(msk), P4a_lag(msk));
disp(['P4:   r = ', num2str(r(1,2)), ', p = ', num2str(p(1,2))])


%% Does it correlate with ox/spice
P4.isoDO(P4.isoDO>700)=NaN;

CS01.oxML(CS01.oxML<10)=CS01.oxML(CS01.oxML<10).*43;
CS01.oxML2(CS01.oxML2<10)=CS01.oxML2(CS01.oxML2<10).*43;
CS01.ox(CS01.ox<10)=CS01.ox(CS01.ox<10).*43;
CS01.ox2(CS01.ox2<10)=CS01.ox2(CS01.ox2<10).*43;

sig_thresh=[26 26.5 26.75];
isoDO=NaN(3,size(CS01.ox,2));
for i =1:length(sig_thresh)
    for ii=1:size(CS01.ox,2)
        ox=mean([CS01.oxML(:,ii) CS01.oxML2(:,ii) CS01.ox(:,ii) CS01.ox2(:,ii)],...
           2,'omitnan');
        if sum(~isnan(ox))>1
            msk=~isnan(CS01.s_t(:,ii)) & ~isnan(ox);
            isoDO(i,ii)=interp1(CS01.s_t(msk,ii),ox(msk),sig_thresh(i));
        end
    end
end

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 14 10],'color','w');

subplot(1,2,1)
hold on 
scatter(alpha_PEW_CS01_26p5,isoDO(2,:),20,'filled','markerfacecolor',rgb_x('rose'));
grid on; axis tight; 
title('CS01');
xlabel('Fraction of PEW');
ylabel('DO (\mumol kg^{-1})');
msk=~isnan(alpha_PEW_CS01_26p5) &~isnan(isoDO(2,:));
coeffs=polyfit(alpha_PEW_CS01_26p5(msk),isoDO(2,msk),1);
x=linspace(min(alpha_PEW_CS01_26p5),max(alpha_PEW_CS01_26p5));
fittedy=polyval(coeffs,x,100);
plot(x,fittedy,'linewidth',2','color',rgb_x('rose'));

subplot(1,2,2)
hold on 
scatter(alpha_PEW_P4_26p5,P4.isoDO(2,:),20,'filled','markerfacecolor',rgb_x('steel blue'));
grid on; axis tight; 
title('P4');
xlabel('Fraction of PEW');
ylabel('DO (\mumol kg^{-1})');
msk=~isnan(alpha_PEW_P4_26p5) &~isnan(P4.isoDO(2,:));
coeffs=polyfit(alpha_PEW_P4_26p5(msk),P4.isoDO(2,msk),1);
x=linspace(min(alpha_PEW_P4_26p5),max(alpha_PEW_P4_26p5));
fittedy=polyval(coeffs,x);
plot(x,fittedy,'linewidth',2','color',rgb_x('steel blue'));

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',10);

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/PEWvsDO.jpg -m4


%% How about hypoxic years

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 10 10],'color','w');
hold on 
scatter(alpha_PEW_CS01_26p5,isoDO(2,:),10,'filled','markerfacecolor',[0.6 0.6 0.6]);
lnes=linspecer(7);
yr=year(CS01.mtime);
count=0;
for i=[2000, 2002, 2005, 2019, 2020, 2021, 2023]
    count=count+1;
    msk=yr==i;
    % if i==2002
    %     keyboard
    % end
    ss(count)=scatter(alpha_PEW_CS01_26p5(msk),isoDO(2,msk),40,'filled', ...
        'markerfacecolor',lnes(count,:),'markeredgecolor','k');
end

grid on; axis tight; 
title('CS01');
xlabel('Fraction of PEW');
ylabel('DO (\mumol kg^{-1})');
legend(ss,{'2000', '2002', '2005', '2019', '2020', '2021', '2023'});

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',10); 

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/CS01_PEWvsDO.jpg -m4

%% Calculate O*

msk = ~isnan(alpha_PEW_CS01_26p5) & ~isnan(isoDO(2,:));
coeffs = polyfit(alpha_PEW_CS01_26p5(msk), isoDO(2,msk), 1);  % [slope, intercept]
O_fit = polyval(coeffs, alpha_PEW_CS01_26p5);                % fitted DO from PEW fraction
CS01.Ostar = isoDO(2,:) - O_fit;                              % residual = observed - fitted

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoDO(2,:));
coeffs = polyfit(alpha_PEW_P4_26p5(msk), P4.isoDO(2,msk), 1);
O_fit = polyval(coeffs, alpha_PEW_P4_26p5);
P4.Ostar = P4.isoDO(2,:) - O_fit;

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 14 10],'color','w');

subplot(1,2,1)
hold on 
scatter(alpha_PEW_CS01_26p5,CS01.Ostar,20,'filled','markerfacecolor',rgb_x('rose'));
grid on; axis tight; 
title('CS01');
xlabel('Fraction of PEW');
ylabel('O^* (\mumol kg^{-1})');
% msk=~isnan(alpha_PEW_CS01_26p5) &~isnan(isoDO(2,:));
% coeffs=polyfit(alpha_PEW_CS01_26p5(msk),isoDO(2,msk),1);
% x=linspace(min(alpha_PEW_CS01_26p5),max(alpha_PEW_CS01_26p5));
% fittedy=polyval(coeffs,x,100);
% plot(x,fittedy,'linewidth',2','color',rgb_x('rose'));

subplot(1,2,2)
hold on 
scatter(alpha_PEW_P4_26p5,P4.Ostar,20,'filled','markerfacecolor',rgb_x('steel blue'));
grid on; axis tight; 
title('P4');
xlabel('Fraction of PEW');
ylabel('O^* (\mumol kg^{-1})');
% msk=~isnan(alpha_PEW_P4_26p5) &~isnan(P4.isoDO(2,:));
% coeffs=polyfit(alpha_PEW_P4_26p5(msk),P4.isoDO(2,msk),1);
% x=linspace(min(alpha_PEW_P4_26p5),max(alpha_PEW_P4_26p5));
% fittedy=polyval(coeffs,x);
% plot(x,fittedy,'linewidth',2','color',rgb_x('steel blue'));

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',10);

%% Plot Ostar time-series

f1 = figure('units','centimeters','outerposition',...
    [0.01 0.01 10 8],'color','w');
hold on;

% --- Raw Ostar points ---
scatter(CS01.mtime, CS01.Ostar, 20, 'filled', '^', ...
    'markerfacecolor', rgb_x('rose'), 'MarkerFaceAlpha', 0.4);
scatter(P4.mtime, P4.Ostar, 20, 'filled', '^', ...
    'markerfacecolor', rgb_x('steel blue'), 'MarkerFaceAlpha', 0.2);

xl=xlim;
line(xl,[0 0],'linestyle','--','color','k')

% --- Annual means: CS01 ---
[Y_cs, ~, ~] = datevec(CS01.mtime);
years_cs = unique(Y_cs);
mean_Ostar_CS01 = nan(size(years_cs));
for k = 1:length(years_cs)
    idx = (Y_cs == years_cs(k));
    mean_Ostar_CS01(k) = nanmean(CS01.Ostar(idx));
end
l2 = plot(datenum(years_cs,1,1), mean_Ostar_CS01, ...
    'color', rgb_x('rose'), 'LineWidth', 1.5);

% --- Annual means: P4 ---
[Y_p4, ~, ~] = datevec(P4.mtime);
years_p4 = unique(Y_p4);
mean_Ostar_P4 = nan(size(years_p4));
for k = 1:length(years_p4)
    idx = (Y_p4 == years_p4(k));
    mean_Ostar_P4(k) = nanmean(P4.Ostar(idx));
end
l1 = plot(datenum(years_p4,1,1), mean_Ostar_P4, ...
    'color', rgb_x('steel blue'), 'LineWidth', 1.5);

% --- Labels and formatting ---
xlabel('Year');
ylabel('O^* (\mumol kg^{-1})');
axdate(20);  % rotate date ticks
box on; grid on; axis tight;
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
legend([l1, l2], {'P4', 'CS01'}, 'Location', 'southwest', 'FontSize', 8);


%% 
export_fig /Users/samuelstevens/Dropbox/Hakai/gyre/supply/figs/OstarTS.jpg -m4

%%

% 
% subplot(2,2,3)
% scatter(alpha_SAW_P4_26p5,P4.isoDO(2,:),20,'filled','markerfacecolor',rgb_x('steel blue'));
% msk=~isnan(P4.isoDO(2,:)) & ~isnan(alpha_SAW_P4_26p5);
% [r,p] = corrcoef(alpha_SAW_P4_26p5(msk),P4.isoDO(2,msk))
% grid on; axis tight; 
% title('P4');

% subplot(1,3,2)
% scatter(alpha_SAW_P4_26p5,P4.isoSpice(2,:),20,'filled','b');
% msk=~isnan(P4.isoSpice(2,:)) & ~isnan(alpha_SAW_P4_26p5);
% [r,p] = corrcoef(alpha_SAW_P4_26p5(msk),P4.isoSpice(2,msk))