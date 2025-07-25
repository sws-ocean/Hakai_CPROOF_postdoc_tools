%% Source water plot
clc
clear

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

C=load('CUCdataR1.mat');
P=load('WPdataR1.mat');
shelf=load('isoDO_shelfMn.mat');

%% Figure

col=[0.9 0.9 0.9];% [255 214 140]/255; % YELLOW!
yrGrid=2003:2023;
for i=1:length(yrGrid)
    decyrs(i)=datenum(yrGrid(i),01,01);
end

odecyrs=[];
oyrGrid=C.cal.calcofiDOa(:,1);
for i=1:length(oyrGrid)
    odecyrs(i)=datenum(oyrGrid(i),01,01);
end

% lnes=linspecer(5);lnes(3,:)=[];

lnes = [ 
    27, 158, 119;
    231, 41, 138;
    217, 95, 2;
    117, 112, 179
] / 255;

%%%%% Figure start
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 9 16],'color',[1 1 1 0.5]);

ax1=axes('position',[0.12 0.55 0.8 0.5]); hold on 
m_proj('equidistant', 'lon', [145 250], 'lat', [24 60]);
[CS,H]=m_etopo2('contourf',[-10000:500:0],'edgecolor','none');
H.FaceAlpha=0.6;
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linest','none','linewidth',0.5,'tickdir','out','xaxislocation','top');
colormap(m_colmap('blues'));

m_scatter(mod(C.cal.allL + 360, 360),C.cal.allLL,2,'filled','markerfacecolor',lnes(1,:),...
    'markeredgecolor',lnes(1,:));
m_scatter(nanmean(mod(C.P4.lon + 360, 360)),nanmean(C.P4.lat),100,'filled','markerfacecolor',...
    lnes(2,:),'marker','p','markeredgecolor','w');
m_scatter(mod(shelf.isoDO.longitude + 360, 360),shelf.isoDO.latitude,...
    4,'filled','markerfacecolor',lnes(3,:));
m_scatter(P.all_l,P.all_ll,2,'filled','markerfacecolor',lnes(4,:));
text(0.025, 0.95, 'a)', 'Units', 'normalized');

%%%%%% Timeseries
ax2=axes('position',[0.12 0.35 0.8 0.3]);  
hold on;
ss={};

ts=smoothdata(C.cal.calcofiDOa(:,2),'movmean',1,'omitnan');
jb(1)=jbfill(odecyrs,ts-C.cal.AseDO'/2,flipud(ts+C.cal.AseDO'/2),...
    lnes(1,:), lnes(1,:), 1, 0.1);
lo=plot(odecyrs, ts, '-', 'color', lnes(1,:),'linewidth',0.5);
scatter(odecyrs,ts, 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(1,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
coeffs=polyfit(odecyrs, C.cal.calcofiDOa(:,2),1);
msk=~isnan(C.cal.calcofiDOa(:,2)) & odecyrs'>=datenum(2003,01,01);
[s,b]=tsreg(odecyrs(msk), C.cal.calcofiDOa(msk,2)');
[slopeBoot, ~] = bootstrapTrend(odecyrs(msk), C.cal.calcofiDOa(msk,2)', 1000);
fittedy=polyval([s b],odecyrs(odecyrs'>=datenum(2003,01,01)));
trndI(4)=plot(odecyrs(odecyrs'>=datenum(2003,01,01)),...
    fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
trnd(4)=plot(odecyrs(odecyrs'>=datenum(2003,01,01)),...
    fittedy,'--','color',lnes(1,:),'linewidth',1);
ss{4}=sprintf('%3.1f\x00B1%2.1f',s(1)*10*365,2*nanstd(slopeBoot)* 10 * 365);

ts=smoothdata(C.P4.isoDOa(1,:),'movmean',1,'omitnan');
jb(2)=jbfill(decyrs,ts-C.P4.AseDO/2,flipud(ts+C.P4.AseDO/2),...
    lnes(2,:), lnes(2,:), 1, 0.1);
lo=plot(decyrs, ts, '-', 'color', lnes(2,:),'linewidth',0.5);
scatter(decyrs,ts, 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(2,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
coeffs=polyfit(decyrs, C.P4.isoDOa(1,:),1);
[s,b]=tsreg(decyrs, C.P4.isoDOa(1,:));
[slopeBoot, ~] = bootstrapTrend(decyrs, shelf.isoDO.AmnDO(2,:), 1000);
fittedy=polyval([s b],decyrs);
trndI(3)=plot(decyrs,fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
trnd(3)=plot(decyrs,fittedy,'--','color',lnes(2,:),'linewidth',1);
ss{3}=sprintf('%3.1f\x00B1%2.1f',s(1)*10*365,2*nanstd(slopeBoot)* 10 * 365);

ts=smoothdata(shelf.isoDO.AmnDO(2,:),'movmean',1,'omitnan');
jb(3)=jbfill(decyrs,ts-shelf.isoDO.AseDO(1,:)/2,flipud(ts+shelf.isoDO.AseDO(1,:)/2),...
    lnes(3,:), lnes(3,:), 1, 0.1);
plot(decyrs,ts, '-', 'color', lnes(3,:),'linewidth',0.5);
scatter(decyrs,ts, 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(3,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
msk=~isnan(shelf.isoDO.AmnDO(2,:));
coeffs=polyfit(decyrs(msk), shelf.isoDO.AmnDO(2,msk),1);
[s,b]=tsreg(decyrs(msk), shelf.isoDO.AmnDO(2,msk));
[slopeBoot, ~] = bootstrapTrend(decyrs(msk), shelf.isoDO.AmnDO(2,msk), 1000);
fittedy=polyval([s b],decyrs);
trndI(2)=plot(decyrs,fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
trnd(2)=plot(decyrs,fittedy,'--','color',lnes(3,:),'linewidth',1);
ss{2}=sprintf('%3.1f\x00B1 5.5',s(1)*10*365);

odecyrs=[];
for i=1:length(P.oyrGrid)
    odecyrs(i)=datenum(P.oyrGrid(i),01,01);
end

ts=smoothdata(P.isoDOa,'movmean',1,'omitnan');
jb(4)=jbfill(odecyrs,inpaint_nans(ts)-inpaint_nans(P.AseDO)/2,...
    flipud(inpaint_nans(ts)+inpaint_nans(P.AseDO)/2),...
    lnes(4,:), lnes(4,:), 1, 0.1);
lo=plot(odecyrs, ts, '-', 'color', lnes(4,:),'linewidth',0.5);
scatter(odecyrs, ts, 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(4,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
coeffs=polyfit(odecyrs(11:end), P.isoDOa(11:end),1);
[s,b]=tsreg(odecyrs(11:end), [P.isoDOa(11:end)]);
fittedy=polyval([s b],odecyrs(11:end));
trndI(1)=plot(odecyrs(11:end),fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
trnd(1)=plot(odecyrs(11:end),fittedy,'--','color',lnes(4,:),'linewidth',1);

x = odecyrs(11:end);
y = P.isoDOa(11:end);
idxs=bootrnd(length(x),10000);
bootSlopes = zeros(1, 10000);
for i = 1:10000
    idx = idxs(:,i);
    bootSlopes(i) = tsreg(x(idx), y(idx));
end
confIntervals = diff(prctile(bootSlopes, [2.5 97.5]))/2*10*365;

ss{1}=sprintf('%3.1f\x00B1%2.1f',s(1)*10*365,confIntervals);

grid on
axis tight;axdate(10);
ylabel('O$_2$ ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
xlabel('Year');
% xticklabels([]);
text(0.025, 0.975, 'b)', 'Units', 'normalized');


%%%%%%%%%%%%%
odecyrs=[];
oyrGrid=C.cal.calcofiDOa(:,1);
for i=1:length(oyrGrid)
    odecyrs(i)=datenum(oyrGrid(i),01,01);
end

% ax2=axes('position',[0.12 0.35 0.8 0.3]);  
ax3=axes('position',[0.12 0.075 0.8 0.2]);  
hold on;

ts=smoothdata(detrend(shelf.isoDO.AmnDO(2,:)),'movmean',3);
scatter(decyrs,detrend(shelf.isoDO.AmnDO(2,:)), 20, 'o', 'filled', 'markerfacecolor', ...
    lnes(3,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.6);
plot(decyrs,ts, '-', 'color', lnes(3,:),'linewidth',1);

odecyrs=[];
for i=1:length(P.oyrGrid)
    odecyrs(i)=datenum(P.oyrGrid(i),01,01);
end

ts=smoothdata(detrend(inpaint_nans(P.isoDOa)),'movmean',3);
scatter(odecyrs+(8*365), detrend(inpaint_nans(P.isoDOa)), 20, 'o', 'filled', 'markerfacecolor', ...
    lnes(4,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.6);
plot(odecyrs(~isnan(P.isoDOa))+(8*365), ts(~isnan(P.isoDOa)), '-', 'color',...
    lnes(4,:),'linewidth',1);

grid on
axis tight;
ylabel('O'' ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
% xlabel('Year');
text(0.025, 0.975, 'c)', 'Units', 'normalized');
xlim([datenum(2001,1,1) datenum(2023,1,1)]);
axdate(4);
xtickangle(45);
xlabel('Year');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);


axes(ax2);
lg=legend(trnd,ss,'fontsize',8);
title(lg,'O$_2$ Trends ($\mu$mol kg$^{-1}$ 10 yrs$^{-1}$)', 'Interpreter',...
    'latex','fontsize',8);
lg.NumColumns = 2;

uistack(jb,'bottom');
uistack(trndI,'top');
uistack(trnd,'top');


set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
axes(ax3)

 %%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/sourceWatersR1.pdf -dpdf -nofontswap
