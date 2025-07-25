%% Climate indices 
clear 
addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

load isoDO_DFOandHakai.mat
load aPDO.mat
load aNPGO


%% 
sig_thresh = [26 26.5 26.75];

isoDO.oxygen(isoDO.oxygen < 15) = isoDO.oxygen(isoDO.oxygen < 15) * 43.7;
isoDO.oxygen(isoDO.oxygen > 600)=NaN;
yrs = 2003:2023;
num_years = length(yrs);
num_sig_thresh = length(sig_thresh);

% Initialize matrices for annual means and 95% CI for depth, temperature, and salinity
isoDO.AmnDO = NaN(num_sig_thresh, num_years);

for i = 2003:2023
    year_idx = i - 2002;
    msk = isoDO.time > datenum(i, 1, 1) & isoDO.time < datenum(i + 1, 1, 1) & isoDO.latitude < 52;

    for j = 1:num_sig_thresh
        % Dissolved Oxygen
        data_DO = isoDO.oxygen(j, msk);
        data_DO = data_DO(~isnan(data_DO)); % Remove NaNs
        if ~isempty(data_DO)
            isoDO.AmnDO(j, year_idx) = mean(data_DO);
        end 
    end
end

aPDO(aPDO(:,1)<2003 | aPDO(:,1)>2023,:)=[];

figure; lnes=linspecer(3);
tmp=detrend(aPDO(:,2));
for i=1:3
    isoDO.AmnDOdt(i,:)=detrend(inpaint_nans(isoDO.AmnDO(i,:))); % DETREND DATA

    subplot(2,3,i)
    scatter(isoDO.AmnDO(i,:),aPDO(:,2),'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:));
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),aPDO(msk',2));
    text(0.9,0.1,[num2str(r(2)),'/',num2str(pval(2))],'units','normalized')
   title('annual shelf vs. PDO')

        subplot(2,3,i+3)
    scatter(isoDO.AmnDOdt(i,:),tmp,'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:));
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),tmp(msk));
    text(0.9,0.1,[num2str(r(2)),'/',num2str(pval(2))],'units','normalized')
    title('DETRENDED: annual shelf vs. PDO')


end
 
% figure;
% for i=1:3
%     subplot(1,3,i)
%     scatter(CS03.isoDOa(i,:),aPDO(pMsk,2),'filled','markerfacecolor',lnes(i,:));
%     msk=~isnan(CS03.isoDOa(i,:));
%     [r,pval]=corrcoef(CS03.isoDOa(i,msk),aPDO(msk));
%     text(0.9,0.1,[num2str(r(2)),'/',num2str(pval(2))],'units','normalized')
%     sgtitle('CS03 vs. PDO')
% end

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 15 12],'color','w');
tmp=detrend(aNPGO);
for i=1:3
    subplot(2,3,i)
    scatter(isoDO.AmnDO(i,:),aNPGO,'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:));
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),aNPGO(msk));
    text(0.99,0.1,{['r=',num2str(r(2),'%2.2f')];...
        ['p-value=',num2str(pval(2),'%2.2f')]},'units','normalized',...
        'horizontalalignment','right');
    title([num2str(sig_thresh(i),'%2.2f') ' \sigma_\theta']);
    xlabel('$\overline{\rm{O}_2}$ ($\mu$mol~kg$^{-1}$)','interpreter','latex');
    ylabel('NPGO Index')
    grid on;

    subplot(2,3,i+3)
    scatter(isoDO.AmnDOdt(i,:),tmp,'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:));
    [r,pval]=corrcoef(isoDO.AmnDOdt(i,msk),tmp(msk));
    text(0.99,0.1,{['r=',num2str(r(2),'%2.2f')];...
        ['p-value=',num2str(pval(2),'%2.2f')]},'units','normalized',...
        'horizontalalignment','right')
    grid on
    xlabel('O'' ($\mu$mol~kg$^{-1}$)','interpreter','latex');
    ylabel('Detrended NPGO Index')
end
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/NPGOvsO.jpg -m4 -nofontswap

%% Correlate climate indices and oxygen at all isopycnals with yearly time lags

% Define time lags
lags = -10:10;
num_lags = length(lags);

% Determine the number of isopycnal levels
num_isopycnals = size(isoDO.oxygen, 1);

% Ensure the climate indices cover the same years
aPDO(aPDO(:,1)<2003 | aPDO(:,1)>2023,:) = [];
PDO = detrend(aPDO(:,2)); % PDO index values for 2003-2023
NPGO = detrend(aNPGO);    % NPGO index values for 2003-2023

% Initialize matrices to store correlation coefficients and p-values
r_PDO = NaN(num_isopycnals, num_lags);
p_PDO = NaN(num_isopycnals, num_lags);
r_NPGO = NaN(num_isopycnals, num_lags);
p_NPGO = NaN(num_isopycnals, num_lags);

% Loop over each time lag
for lag_idx = 1:num_lags
    lag = lags(lag_idx);

    % Adjust the indices based on the lag
    if lag > 0
        % Oxygen lags climate index
        oxygen_ts = isoDO.AmnDOdt(:, 1:end - lag);
        idx_PDO_lagged = PDO((1+lag):end);
        idx_NPGO_lagged = NPGO((1+lag):end);
    elseif lag < 0
        % Climate index lags oxygen
        oxygen_ts = isoDO.AmnDOdt(:, (1 - lag):end);
        idx_PDO_lagged = PDO(1:end + lag);
        idx_NPGO_lagged = NPGO(1:end + lag);
    else
        % Zero lag
        oxygen_ts = isoDO.AmnDOdt;
        idx_PDO_lagged = PDO;
        idx_NPGO_lagged = NPGO;
    end

    % Number of data points after lag adjustment
    num_pts = size(oxygen_ts, 2);

    % Compute correlations for each isopycnal level
    for j = 1:num_isopycnals
        oxy = oxygen_ts(j, :)';
        % Correlation with PDO
        msk_PDO = ~isnan(oxy) & ~isnan(idx_PDO_lagged);
        if sum(msk_PDO) >= 2
            [r, pval] = corrcoef(oxy(msk_PDO), idx_PDO_lagged(msk_PDO));
            r_PDO(j, lag_idx) = r(1, 2);
            p_PDO(j, lag_idx) = pval(1, 2);
        end
        % Correlation with NPGO
        msk_NPGO = ~isnan(oxy) & ~isnan(idx_NPGO_lagged);
        if sum(msk_NPGO) >= 2
            [r, pval] = corrcoef(oxy(msk_NPGO), idx_NPGO_lagged(msk_NPGO));
            r_NPGO(j, lag_idx) = r(1, 2);
            p_NPGO(j, lag_idx) = pval(1, 2);
        end
    end
end

% Plot the correlation coefficients as a function of isopycnal levels and time lags
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 20 12],'color','w');
subplot(1,2,1);
imagesc(lags, 1:num_isopycnals, p_PDO);
colorbar;
xlabel('Time Lag (years)');
ylabel('\sigma_\theta (kg m^{-3})');
yticks([1,2,3]);
yticklabels({'26.00','26.50','26.75'});
title('P values between O'' and detrended PDO', ...
    'interpreter','latex');

subplot(1,2,2);
imagesc(lags, 1:num_isopycnals, r_PDO);
colorbar;
xlabel('Time Lag (years)');
ylabel('\sigma_\theta (kg m^{-3})');
yticks([1,2,3]);
yticklabels({'26.00','26.50','26.75'});
title('R values between O'' and detrended PDO', ...
    'interpreter','latex');
% colormap(m_colmap(('jet')));
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/PDOvsO.jpg -m4 -nofontswap

% Example for NPGO:
figure;
subplot(1,2,1)
imagesc(lags, 1:num_isopycnals, p_NPGO);
colorbar;
xlabel('Time Lag (years)');
ylabel('Isopycnal Level Index');
title('P values between Oxygen and NPGO');

subplot(1,2,2)
imagesc(lags, 1:num_isopycnals, r_NPGO);
colorbar;
xlabel('Time Lag (years)');
ylabel('Isopycnal Level Index');
title('Correlation Coefficient between Oxygen and NPGO');


%% Integrated PDO


iPDO=cumsum(aPDO(:,2));
figure;
plot(2003:2023,smoothdata(isoDO.AmnDO,2,'movmean',5,'omitnan')-nanmean(isoDO.AmnDO,2));
yyaxis right
plot(2003:2023,iPDO); 

