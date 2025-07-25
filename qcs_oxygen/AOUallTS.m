%% load CS09 and P4 timeseries and look at trend
clc
clear

load isoDO_shelfMn.mat
WP=load('WPdataR1.mat');

sig_thresh = [26 26.5 26.7];
yrGrid=2003:2023;
load CS09_ctd.mat
CS09.ox(:,[1:4 10])=[];
CS09.s_t(:,[1:4 10])=[];
CS09.mtime([1:4 10])=[];
CS09.aou(:,[1:4 10])=[];


CS09.isoDO=NaN(3,length(CS09.mtime));
CS09.isoAOU=NaN(3,length(CS09.mtime));
CS09.isoZ=NaN(3,length(CS09.mtime));
for i=1:length(sig_thresh)
    for ii=1:length(CS09.isoDO)
        msk=~isnan(CS09.s_t(:,ii));
        if sum(msk)>0
            CS09.isoDO(i,ii)=interp1(CS09.s_t(msk,ii),CS09.ox(msk,ii),sig_thresh(i));
            CS09.isoAOU(i,ii)=interp1(CS09.s_t(msk,ii),CS09.aou(msk,ii),sig_thresh(i));
            
            if CS09.isoAOU(i,ii)==0 && i==2
                keyboard
            end

            CS09.isoZ(i,ii)=interp1(CS09.s_t(msk,ii),CS09.pr(msk,ii),sig_thresh(i));
        end


    end
end

CS09.isoDOa=NaN(3,length(yrGrid));
CS09.isoAOUa=NaN(3,length(yrGrid));
CS09.AseDO=NaN(3,length(yrGrid));
CS09.deepST=NaN(3,length(yrGrid));
CS09.deepO=NaN(3,length(yrGrid));
CS09.deepOe=NaN(3,length(yrGrid));

yrs=year(CS09.mtime);
mnth=month(CS09.mtime);
for i=yrGrid(1):yrGrid(end)
    msk=yrs==i & mnth>=8 & mnth<=10;
    CS09.isoDOs(:,i-2002)=mean(CS09.isoDO(:,msk),2,'omitnan');
    CS09.deepST(i-2002)=mean(CS09.s_t(180,msk),'omitnan');
    CS09.deepO(i-2002)=mean(CS09.ox(180,msk),'omitnan');
    CS09.deepOe(:,i-2002)=mean(abs(CS09.ox(180,msk)-CS09.deepO(i-2002)),2,'omitnan');

    msk=yrs==i;
    CS09.isoDOa(:,i-2002)=mean(CS09.isoDO(:,msk),2,'omitnan');
    % CS09.AseDO(:,i-2002)=mean(abs(CS09.isoDO(:,msk)-CS09.isoDOa(:,i-2002)),2,'omitnan');
    CS09.AseDO(:,i-2002)=2*std(CS09.isoDO(:,msk),0,2,'omitnan');
    CS09.isoAOUa(:,i-2002)=mean(CS09.isoAOU(:,msk),2,'omitnan');
    CS09.AOUe(:,i-2002)=2*std(CS09.isoAOU(:,msk),0,2,'omitnan');
end


% CS01
load CS01_ctd.mat
msk=CS01.lat>51;
CS01.lon(msk)=NaN;CS01.lat(msk)=NaN; 
CS01.ox(msk)=NaN;

CS01.isoDO=NaN(3,length(CS01.mtime));
CS01.isoAOU=NaN(3,length(CS01.mtime));
CS01.isoZ=NaN(3,length(CS01.mtime));

for i=1:length(sig_thresh)
    for ii=1:length(CS01.isoDO)
        msk=~isnan(CS01.s_t(:,ii));
        if sum(msk)>0
            CS01.isoDO(i,ii)=interp1(CS01.s_t(msk,ii),CS01.ox(msk,ii),sig_thresh(i));
            CS01.isoAOU(i,ii)=interp1(CS01.s_t(msk,ii),CS01.aou(msk,ii),sig_thresh(i));
            CS01.isoZ(i,ii)=interp1(CS01.s_t(msk,ii),CS01.pr(msk,ii),sig_thresh(i));
        end
    end
end

CS01.isoDOa=NaN(3,length(yrGrid));
CS01.isoAOUa=NaN(3,length(yrGrid));
CS01.AseDO=NaN(3,length(yrGrid));
CS01.deepST=NaN(3,length(yrGrid));
CS01.deepO=NaN(3,length(yrGrid));
CS01.deepOe=NaN(3,length(yrGrid));

yrs=year(CS01.mtime);
mnth=month(CS01.mtime);
for i=yrGrid(1):yrGrid(end)
    msk=yrs==i & mnth>=9 & mnth<=9;
    CS01.isoDOs(:,i-2002)=mean(CS01.isoDO(:,msk),2,'omitnan');
    CS01.deepST(i-2002)=mean(CS01.s_t(180,msk),'omitnan');
    CS01.deepO(i-2002)=mean(CS01.ox(180,msk),'omitnan');
    CS01.deepOe(:,i-2002)=mean(abs(CS01.ox(180,msk)-CS01.deepO(i-2002)),2,'omitnan');

    msk=yrs==i;
    CS01.isoDOa(:,i-2002)=mean(CS01.isoDO(:,msk),2,'omitnan');
    % CS01.AseDO(:,i-2002)=mean(abs(CS01.isoDO(:,msk)-CS01.isoDOa(:,i-2002)),2,'omitnan');
    CS01.AseDO(:,i-2002)=2*std(CS01.isoDO(:,msk),0,2,'omitnan');

    CS01.isoAOUa(:,i-2002)=mean(CS01.isoAOU(:,msk),2,'omitnan');
    CS01.AOUe(:,i-2002)=2*std(CS01.isoAOU(:,msk),0,2,'omitnan');
end


%% Plot

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 11 11],'color',[1 1 1]);

subplot(5,1,2:5)

lnes=[1 0 0;0 0 0;rgb_x('royal purple')];

jb(1)=jbfill(yrGrid,inpaint_nans(CS01.isoAOUa(2,:)'-CS01.AOUe(2,:)'),...
    flipud(inpaint_nans(CS01.isoAOUa(2,:)'+CS01.AOUe(2,:)')),...
    lnes(1,:), lnes(1,:), 1, 0.1);
lo=plot(yrGrid,CS01.isoAOUa(2,:), '-', 'color', lnes(1,:),'linewidth',0.5);
trnd(3)=scatter(yrGrid,CS01.isoAOUa(2,:), 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(1,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
msk=~isnan(CS01.isoAOUa(2, :));
[s,b]=tsreg(yrGrid(msk), CS01.isoAOUa(2,msk));
[slopeBoot, ~] = bootstrapTrend(yrGrid(msk), CS01.isoAOUa(2,msk), 1000);
fittedy=polyval([s b],yrGrid);
trndI(3)=plot(yrGrid,fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
plot(yrGrid,fittedy,'--','color',lnes(1,:),'linewidth',1);
ss{3}=sprintf('CS01: %3.1f\x00B1%2.1f',s(1)*10,2*nanstd(slopeBoot)* 10);

jb(1)=jbfill(yrGrid,inpaint_nans(CS09.isoAOUa(2,:)'-CS09.AOUe(2,:)'),...
    flipud(inpaint_nans(CS09.isoAOUa(2,:)'+CS09.AOUe(2,:)')),...
    lnes(2,:), lnes(2,:), 1, 0.1);
lo=plot(yrGrid,CS09.isoAOUa(2,:), '-', 'color', lnes(2,:),'linewidth',0.5);
trnd(2)=scatter(yrGrid,CS09.isoAOUa(2,:), 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(2,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
msk=~isnan(CS09.isoAOUa(2, :));
[s,b]=tsreg(yrGrid(msk), inpaint_nans(CS09.isoAOUa(2,msk)));
[slopeBoot, ~] = bootstrapTrend(yrGrid(msk), CS09.isoAOUa(2,msk), 1000);
fittedy=polyval([s b],yrGrid);
trndI(2)=plot(yrGrid,fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
plot(yrGrid,fittedy,'--','color',lnes(2,:),'linewidth',1);
ss{2}=sprintf('CS09: %3.1f\x00B1%2.1f',s(1)*10,2*nanstd(slopeBoot)* 10);

jb(1)=jbfill(WP.oyrGrid,inpaint_nans(WP.isoAOUa'-WP.AseAOU'),...
    flipud(inpaint_nans(WP.isoAOUa'+WP.AseAOU')),...
    lnes(3,:), lnes(3,:), 1, 0.1);
lo=plot(WP.oyrGrid,WP.isoAOUa, '-', 'color', lnes(3,:),'linewidth',0.5);
trnd(1)=scatter(WP.oyrGrid,WP.isoAOUa, 10, 'o', 'filled', 'markerfacecolor', ...
    lnes(3,:), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.9);
msk=WP.oyrGrid>=2003;
[s,b]=tsreg(WP.oyrGrid(msk), inpaint_nans(WP.isoAOUa(msk)));
[slopeBoot, ~] = bootstrapTrend(WP.oyrGrid(msk), WP.isoAOUa(msk), 1000);
fittedy=polyval([s b],yrGrid);
trndI(1)=plot(yrGrid,fittedy,'-','color',[1 1 1 0.5],'linewidth',2);
plot(yrGrid,fittedy,'--','color',lnes(3,:),'linewidth',1);
ss{1}=sprintf('NWP: %3.1f\x00B1%2.1f',s(1)*10,2*nanstd(slopeBoot)* 10);

lg=legend(trnd,ss,'fontsize',7);
title(lg,'AOU Trends ($\mu$mol kg$^{-1}$ 10 yrs$^{-1}$)', 'Interpreter',...
    'latex','fontsize',10);
lg.NumColumns = 2;

axis tight;
grid on;
ylabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');

uistack(jb,'bottom');
uistack(trndI,'top');
uistack(trnd,'top');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');

%%

export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/AOUtsR1.pdf -dpdf -nofontswap
