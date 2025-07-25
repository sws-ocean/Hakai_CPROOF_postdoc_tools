%% AOU on shelf
addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear
load Dcorrected_allProfs_20241011.mat
load glider_allProfs.mat canyonX canyonZ canyonAX
load S2.mat 
load H1.mat
M.Scott2=S2;  M.Hak1=H1;
clearvars H1 S2
% load DO_1990-2024_QCS.mat
% load DO_1990-2024_QCS_v2.mat
load isoDO_DFOandHakai.mat
load CS09_ctd.mat

bottomDO.oxygen(bottomDO.oxygen>600)=NaN;

sig_thresh = [26 26.5 26.75];

allProfs.s_t=gsw_sigma0(allProfs.SA,allProfs.CT);
allProfs.spice=gsw_spiciness0(allProfs.SA,allProfs.CT);


%% Look at along isopycnal AOU gradient

allProfs.AOU=aou(allProfs.salinity,allProfs.potential_temperature,...
    allProfs.oxygen_corrected);
allProfs.isoAOU=NaN(3,length(allProfs.AOU));
allProfs.isoDO=NaN(3,length(allProfs.AOU));
allProfs.isoSpice=NaN(3,length(allProfs.AOU));

sig_thresh=[26 26.5 26.7];

for i=1:size(allProfs.AOU,2)
        for ii=1:length(sig_thresh)
            msk=~isnan(allProfs.AOU(:,i));
            if sum(msk)>1
                allProfs.isoAOU(ii,i)=interp1(inpaint_nans(allProfs.s_t(msk,i)),allProfs.AOU(msk,i),...
                    sig_thresh(ii));
                allProfs.isoDO(ii,i)=interp1(inpaint_nans(allProfs.s_t(msk,i)),allProfs.oxygen_corrected(msk,i),...
                    sig_thresh(ii));
                allProfs.isoSpice(ii,i)=interp1(inpaint_nans(allProfs.s_t(msk,i)),allProfs.spice(msk,i),...
                    sig_thresh(ii));
            end
        end
end


%% make along canyon average

xGrid=0.25:0.25:159.75;
isoAOUmn=NaN(3,length(xGrid));
isoAOUsem=NaN(3,length(xGrid));
ctdAOUmn=NaN(3,length(xGrid));
ctdAOUsem=NaN(3,length(xGrid));

load isoDO_DFOandHakai.mat
isoDO.AOU(isoDO.AOU>275)=NaN;

count=0;
for i=0:0.25:159.5
    count=count+1;
    msk=allProfs.canyonX>=i & allProfs.canyonX<=i+0.25;
    
    isoAOUmn(:,count)=mean(allProfs.isoAOU(:,msk),2,'omitnan');
    isoAOUsem(:,count)=2.*std(allProfs.isoAOU(:,msk),0,2,'omitnan');%./...
        % sqrt(sum(~isnan(allProfs.isoAOU(:,msk)),2,'omitnan'));

    msk=isoDO.OSdist/1e3>=i & isoDO.OSdist/1e3<=i+0.25;
    ctdAOUmn(:,count)=mean(isoDO.AOU(:,msk),2,'omitnan');
    ctdAOUsem(:,count)=2.*std(isoDO.AOU(:,msk),0,2,'omitnan');%./...
        % sqrt(sum(~isnan(isoDO.AOU(:,msk)),2,'omitnan'));

end
%%
xGridS=(xGrid-120).*-1;
canyonX=(canyonX-120).*-1;

figure('units','centimeters','outerposition',[0 0 12 15],'color','w');
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
load('BCcoast');
lat_lim=[50.5 52.5]; lon_lim=[-130.5 -127.5];

fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

col=[255 214 140]/255; % YELLOW!
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
deep=-4000:1000:-1000; shelf=[-600:50:-100 0];
[CS,CH]=m_contourf(lon(ilon),lat(ilat),Z',[deep shelf],'linestyle','none');

g=cmocean('gray',600);g(1:390,:)=[];g(end-20:end,:)=[];
cm=[g;m_colmap('blues',50)];
colormap(cm);

ax = colorbar('Position', [0.9 0.8 0.01 0.1]); 
set(ax, 'ylim', [-375 0],'box','off', 'YColor', 'k');
ax.Ticks=[-300  -200  -100     0];
ax.TickLabels={'300';'200';'100';'0'};
title(ax,'$h$ (m)','Interpreter','latex')
CH.FaceAlpha=0.9;

m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',0.5,'tickdir','out','xtick',3,'clipping','off');

tmpAx=((isoDO.OSdist/1e3)-120).*-1;
msk=find(tmpAx>-30 & tmpAx<120);
cmap = cmocean('speed', 256);
cmap(1:30,:)=[];
tmpAx_normalized = (tmpAx(msk) - min(tmpAx(msk))) / (max(tmpAx(msk)) - min(tmpAx(msk)));
colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, tmpAx_normalized);


% Loop through each point and plot it individually
for i = 1:length(msk)
    m_scatter(isoDO.longitude(msk(i)), isoDO.latitude(msk(i)), ...
        15, '^', 'filled', 'MarkerFaceColor', colors(i, :),...
        'markerfacealpha',0.7);
end
set(gca,'clipping','off')

m_plot(canyonAX(:,1),canyonAX(:,2),'k','linewidth',1);

% Build axis
lnes=linspecer(length(0:25:max(canyonX)));
count=0;

for i=-30:30:max(canyonX)
    [~,idx]=min(abs(canyonX-i));
    count=count+1;
    m_scatter(canyonAX(idx,1),canyonAX(idx,2),30,'filled','markerfacecolor','k',...
        'markeredgecolor','w','marker','s');
    % m_text(canyonAX(idx,1)-0.4,canyonAX(idx,2),sprintf('%3.0f km',canyonX(idx)));

end
text(0.05, 0.9, 'a)', 'Units', 'normalized');

nexttile;
hold on;
lnes = linspecer(5);lnes(3:4,:)=[]; lnes(1,:)=[0.5 0.5 0.5]; lnes=flipud(lnes);
coeffs=[]; 
coeffs23=[];
coeffsC=[];
fittedyC=[];
for i=1:length(sig_thresh)
    fittedy=[];
    msk=~isnan(isoAOUmn(i,:));
    % errorbar(xGridS,ctdAOUmn(i,:), ctdAOUsem(i,:), 'LineWidth', 0.25, ...
    %     'CapSize', 2, 'Color', [lnes(i,:)-[0.2 0.2 0.2] 0.025], 'LineStyle', 'none');
    fill([xGridS(msk) fliplr(xGridS(msk))]',...
        [isoAOUmn(i,msk)-isoAOUsem(i,msk) fliplr(isoAOUmn(i,msk)+isoAOUsem(i,msk))]',...
        lnes(i,:),'edgecolor','none','facealpha',0.1);

    scatter(xGridS,isoAOUmn(i,:),7,'filled','markerfacecolor',lnes(i,:),...
        'markeredgecolor','none','markerfacealpha',0.3);
    
    scatter(xGridS,ctdAOUmn(i,:),10,'filled','markerfacecolor',lnes(i,:)-[0.1 0.1 0.1],...
        'markeredgecolor','none','marker','^','markerfacealpha',0.7);

    trendMsk=~isnan(isoAOUmn(i,:)) & xGrid<=120;
    [coeffs(i,1),coeffs(i,2)]=tsreg(xGridS(trendMsk), isoAOUmn(i,trendMsk));
    fittedy(i,:)=polyval(coeffs(i,:),xGridS(trendMsk));
        
    plot(xGridS(trendMsk),fittedy(i,:),'color','w','LineWidth',3);
    sp(i)=plot(xGridS(trendMsk),fittedy(i,:),'color',lnes(i,:),'LineWidth',1.5);

    % trendMsk=~isnan(ctdAOUmn(i,:)) & xGrid<=120;
    % [coeffsC(i,1),coeffsC(i,2)]=tsreg(xGridS(trendMsk), ctdAOUmn(i,trendMsk));
    % fittedyC(i,:)=polyval(coeffsC(i,:),xGridS(xGrid<=120));
    % 
    % plot(xGridS(xGrid<=120),fittedyC(i,:),'color','w','LineWidth',3);
    % plot(xGridS(xGrid<=120),fittedyC(i,:),'--','color',lnes(i,:),'LineWidth',1.5);
    % 
end

%set(gca,'xdir','reverse');
xlabel('Distance from 300 m isobath (km)');
ylabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
grid on; axis tight; ylim([80 270]);
text(0.05, 0.9, 'b)', 'Units', 'normalized');
xticks(-30:30:120);

for i=1:3
    s(i)=scatter(NaN,NaN,10,'filled','markerfacecolor',lnes(i,:),...
        'markeredgecolor','none');
end

legend(s,'26.00 \sigma_\theta', '26.50 \sigma_\theta', '26.75 \sigma_\theta',...
    'edgecolor','none','color','none')


axC=axes('position',[0.725 0.575 0.1 0.175]); hold on
set(axC, 'Visible', 'off'); % Turn off the entire axes

for i=1:size(cmap,1)
    scatter(1,1-i*0.01,15,'s','markerfacecolor',cmap(i,:),'markeredgecolor','none');
    if i==1
        text(1.2,1-i*0.01,'0','fontsize',8);
    elseif i==round(length(cmap)/2)
        text(1.2,1-i*0.01,'60','fontsize',8);
    elseif i==length(cmap)
        text(1.2,1-i*0.01,'120','fontsize',8);
    end
end
text(1,1+0.01*10,{'Distance from';'300 m isobath (km)'},'fontsize',8,...
    'horizontalalignment','left');


set(findall(gcf,'-property','fontsize'),'fontsize',10);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

% Rough transit time estimate
% 
load SWSadcpS2.mat

msk1=month(adcp{1, 1}.Time)>=5 & month(adcp{1, 1}.Time)<=10;
msk2=adcp{1, 1}.bin_depth>100;
meanAcrossU=mean(adcp{1,1}.velcs(msk1,msk2),2,'omitnan');
% TT(1)=(120e3/max(meanAcrossU))/86400; %days
TT=(120e3/mean(meanAcrossU,'omitnan'))/86400; %days
std(meanAcrossU,'omitnan')
% TT(2)=(120e3/min(meanAcrossU))/86400; %days

OUR=(coeffs(1,1).*120)./TT
OUR=(coeffs(3,1).*120)./TT

%%
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUvsXv6.pdf', ...
%     'ContentType', 'vector', 'Resolution', 100);
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUvsX_R1.pdf -dpdf -nofontswap
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/poster/AOUshelf.png -m6 -nofontswap -transparent

%% Trend uncertainty
num_bootstraps = 1000;
BS = zeros(3,1); % Preallocate for bootstrap uncertainty

for i = 1:3
    trendMsk = ~isnan(isoAOUmn(i,:)) & xGrid <= 120;
    x = xGridS(trendMsk);
    y = isoAOUmn(i,trendMsk);

    coeffsBS = zeros(num_bootstraps, 2); % Store slope and intercept

    for ii = 1:num_bootstraps
        idx = randi(length(x), length(x), 1);
        coeffsBS(ii,:) = tsreg(x(idx)', y(idx)'); % Ensure both slope & intercept are stored
    end

    % Uncertainty in slope (2*std for approx. 95% CI)
    CI = prctile(coeffsBS(:,1), [2.5 97.5]);  % 95% confidence interval for the slope
    uncertainty = diff(CI) / 2;  % Half-width of CI
    BS(i) = uncertainty;
end
