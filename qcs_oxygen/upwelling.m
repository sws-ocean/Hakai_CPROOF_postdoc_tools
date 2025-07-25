%% Bakun
addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear


load S2.mat
load H1.mat
M=struct();
M.S2=S2; M.H1=H1;
clearvars -except M
load S2_17.mat
M.S2_17=S2;
clearvars -except M


M.H1.CTD133m.sigma_theta(32093:35236)=NaN;

%% pull in Bakun
load BAKmay2024.mat
BAK.B(end)=[]; BAK.mtime(end)=[];
dday=day(datetime(datestr(BAK.mtime)),'dayofyear');

BAK.Cmn=NaN(1,365);
BAK.Csig=NaN(1,365);

for i=1:365
    msk=dday==i;
    BAK.Cmn(i)=mean(BAK.B(msk),'omitnan');
    BAK.Csig(i)=2.*std(BAK.B(msk),'omitnan');
end

BAK.B15=BAK.B(year(BAK.mtime)==2015);
BAK.B22=BAK.B(year(BAK.mtime)==2022);
BAK.B23=BAK.B(year(BAK.mtime)==2023);

yrs=year(BAK.mtime);
CUI=NaN(length(unique(yrs))-1,365);
CUIu=NaN(length(unique(yrs))-1,365);
CUIw=NaN(length(unique(yrs))-1,1);
lns=linspecer(57);
figure; hold on;

for i=1967:2023
    for ii=1:365
        msk=yrs==i & dday==ii;
        CUI(i-1966,ii)=mean(BAK.B(msk),'omitnan');
    end
    tmp=CUI(i-1966,:);tmp(tmp<0)=0;
    CUIu(i-1966,:)=cumsum(tmp,'omitnan');
    CUIw(i-1966)=sum(tmp(1:90),'omitnan');

    CUI(i-1966,:)=cumsum(CUI(i-1966,:),'omitnan');
    plot(1:365,CUI(i-1966,:),'color',[lns(i-1966,:) 0.3]);
end

l(1)=plot(1:365,cumsum(BAK.B15,'omitnan'),'g','linewidth',2);
l(2)=plot(1:365,cumsum(BAK.B22,'omitnan'),'r','linewidth',2);
l(3)=plot(1:365,cumsum(BAK.B23,'omitnan'),'b','linewidth',2);
l(4)=plot(mean(CUI,'omitnan'),'k','LineWidth',2);

legend(l,'2015','-400','2023','Mean')
ylabel({'CUI';'(m^3 s^{-1} 100m^{-1})'});
xlabel ('Julian day');

% for i=1:size(CUIu,1)
%     UPd(i)=find(CUIu(i,:)>2000,1,'first');
% end
% 
% figure;
% scatter(1967:2023,UPd,20,'x')
% hold on
% plot(1967:2023,smoothdata(UPd,'movmean',10));


%% Find extreme events in mooring

syr=year(M.S2.time); smnth=month(M.S2.time);
hyr=year(M.H1.time); hmnth=month(M.H1.time);
ddayH=day(datetime(datestr(M.H1.time)),'dayofyear');
ddayS=day(datetime(datestr(M.S2.time)),'dayofyear');

pr=ones(size(M.S2.CTD280m.salinity)).*(280*1.0051);
SA=gsw_SA_from_SP(M.S2.CTD280m.salinity,...
    pr,M.S2.longitude,M.S2.latitude);
PT=gsw_pt0_from_t(SA,M.S2.CTD280m.temperature,pr);
M.S2.CTD280m.AOU=aou(M.S2.CTD280m.salinity,PT,M.S2.CTD280m.oxygen);

pr=ones(size(M.H1.CTD133m.salinity)).*(133*1.0051);
SA=gsw_SA_from_SP(M.H1.CTD133m.salinity,...
    pr,M.H1.longitude,M.H1.latitude);
PT=gsw_pt0_from_t(SA,M.H1.CTD133m.temperature,pr);
M.H1.CTD133m.AOU=aou(M.H1.CTD133m.salinity,PT,M.H1.CTD133m.oxygen);

SclimWM=NaN(1,365);
HclimWM=NaN(1,365);
SclimWA=NaN(1,365);
HclimWA=NaN(1,365);
HrC=NaN(1,365);
SrC=NaN(1,365);
SclimWE=NaN(1,365);
HclimWE=NaN(1,365);
HrCE=NaN(1,365);
SrCE=NaN(1,365);
SclimAE=NaN(1,365);
HclimAE=NaN(1,365);

for i=1:365
    SclimWM(i)=mean(M.S2.CTD280m.oxygen(ddayS==i),'omitnan');
    HclimWM(i)=mean(M.H1.CTD133m.oxygen(ddayH==i),'omitnan');

    SclimWE(i)=std(M.S2.CTD280m.oxygen(ddayS==i),0,'omitnan');
    HclimWE(i)=std(M.H1.CTD133m.oxygen(ddayH==i),0,'omitnan');

    HrC(i)=mean(M.H1.CTD133m.sigma_theta(ddayH==i),'omitnan');
    SrC(i)=mean(M.S2.CTD280m.sigma_theta(ddayH==i),'omitnan');

    HrCE(i)=std(M.H1.CTD133m.sigma_theta(ddayH==i),0,'omitnan');
    SrCE(i)=std(M.S2.CTD280m.sigma_theta(ddayH==i),0,'omitnan');

    SclimWA(i)=mean(M.S2.CTD280m.AOU(ddayS==i),'omitnan');
    HclimWA(i)=mean(M.H1.CTD133m.AOU(ddayH==i),'omitnan');
    
    SclimAE(i)=std(M.S2.CTD280m.AOU(ddayS==i),0,'omitnan');
    HclimAE(i)=std(M.H1.CTD133m.AOU(ddayH==i),0,'omitnan');
end

% Create dummy timeseries
Hdummy=[repmat(HclimWM,1,7);repmat(HclimWM-1*HclimWE,1,7);repmat(HclimWM+1*HclimWE,1,7);...
    repmat(HclimWM-2*HclimWE,1,7);repmat(HclimWM+2*HclimWE,1,7)];
Sdummy=[repmat(SclimWM,1,7);repmat(SclimWM-1*SclimWE,1,7);repmat(SclimWM+1*SclimWE,1,7);
    repmat(SclimWM-2*SclimWE,1,7);repmat(SclimWM+2*SclimWE,1,7)];

HSTdummy=[repmat(HrCE,1,7);repmat(HrC-1*HrCE,1,7);repmat(HrC+1*HrCE,1,7);...
    repmat(HrC-2*HrCE,1,7);repmat(HrC+2*HrCE,1,7)];
SSTdummy=[repmat(SrC,1,7);repmat(SrC-1*SrCE,1,7);repmat(SrC+1*SrCE,1,7);
    repmat(SrC-2*SrCE,1,7);repmat(SrC+2*SrCE,1,7)];

dummyT=datenum(2017,01,01):datenum(2024,01,01)-2;

tmp=smoothdata(M.H1.CTD133m.oxygen,'movmean',24*7,'omitnan');
msk2=find(tmp<interp1(dummyT,Hdummy(4,:),M.H1.time) & tmp<61);

Hso=M.H1.CTD133m.oxygen(msk2);
Hst=M.H1.CTD133m.sigma_theta(msk2);
Ha=M.H1.CTD133m.AOU(msk2);
Ht=M.H1.CTD133m.temperature(msk2);
Hs=M.H1.CTD133m.salinity(msk2);

msk=M.S2.time>datenum(2023,05,01) & M.S2.time<datenum(2023,08,01) & ...
    ~isnan(M.S2.CTD280m.oxygen) & ~isnan(M.S2.CTD280m.sigma_theta);
Sso=M.S2.CTD280m.oxygen(msk); Sst=M.S2.CTD280m.sigma_theta(msk);
Sao=M.S2.CTD280m.aou(msk);
k = boundary(double(Sso'),double(Sst'));

%% interpolate grids

% daily mean
MgridS=ceil(M.S2.time(1)):floor((M.S2.time(end)));
S2I.dMnSt=NaN(length(MgridS),1);

for i=1:length(Mgrid)
    msk=floor(M.S2.time)==MgridS(i);
    S2I.dMnSt(i)=mean(M.S2.CTD280m.sigma_theta(msk));
end

% daily mean
MgridH=ceil(M.H1.time(1)):floor((M.H1.time(end)));
H1I.dMnSt=NaN(length(MgridH),1);

for i=1:length(Mgrid)
    msk=floor(M.H1.time)==MgridH(i);
    H1I.dMnSt(i)=mean(M.H1.CTD133m.sigma_theta(msk));
end


S2I.s_t=interp1(MgridS,S2I.dMnSt,BAK.mtime);
H1I.s_t=interp1(MgridH,H1I.dMnSt,BAK.mtime);


[N, edgesX, edgesY] = histcounts2(S2I.s_t,BAK.B);
figure
set(gca,'TickDir','out'); 
hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.6;shading flat;




%% Make Bakun plot

lnes=[rgb_x('steel blue');rgb_x('rose')];
figure('units','centimeters','outerposition',[0 0 11 16],'color','w');

t = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
axB22=nexttile;
set(gca,'TickDir','out'); 
hold on;

fill([1:365 fliplr(1:365)]',[BAK.Cmn-BAK.Csig fliplr(BAK.Cmn+BAK.Csig)]',...
    [0.9 0.9 0.9],'edgecolor','none');
line([1 365],[0 0],'linestyle',':','color','k','linewidth',0.5);
plot(1:365,smoothdata(BAK.Cmn,'movmean',14),'k--','linewidth',1);

winB=smoothdata(BAK.B22,1,'movmean',14,'omitnan');
tmp=winB;tmp(tmp>0)=0;
plot(1:365,tmp,'b-','LineWidth',0.5);
tmp=winB;tmp(tmp<0)=0;
plot(1:365,tmp,'r-','LineWidth',0.5);
xl=xlim; line(xl,[0,0],'color',[0.8 0.8 0.8],'linewidth',0.75);
grid on;
yticks(-400:100:300); 
ylim([-200 200]);
axis tight; 
xlim([0 365]);
xticks([0 182 365]);
%set(gca,'XAxisLocation','top');
xlabel('Julian day');

ylabel({'Bakun index';'(m^3 s^{-1} 100m^{-1})'});

% text([121 152 182 213 244 274],[-155 -155 -155 -155 -155 -155],...
%     {'M';'J';'J';'A';'S';'O'},'fontsize',6,...
% %     'HorizontalAlignment','center');


axB23=nexttile;
set(gca,'TickDir','out'); 
hold on;

fill([1:365 fliplr(1:365)]',[BAK.Cmn-BAK.Csig fliplr(BAK.Cmn+BAK.Csig)]',...
    [0.9 0.9 0.9],'edgecolor','none');
line([1 365],[0 0],'linestyle',':','color','k','linewidth',0.5);
plot(1:365,smoothdata(BAK.Cmn,'movmean',14),'k--','linewidth',1);

winB=smoothdata(BAK.B23,1,'movmean',14,'omitnan');
tmp=winB;tmp(tmp>0)=0;
plot(1:365,tmp,'b-','LineWidth',0.5);
tmp=winB;tmp(tmp<0)=0;
plot(1:365,tmp,'r-','LineWidth',0.5);
xl=xlim; line(xl,[0,0],'color',[0.8 0.8 0.8],'linewidth',0.75);
grid on;
yticks(-400:100:300); 
ylim([-200 200]);
axis tight;
xlim([0 365]);
xticks([0 182 365]);
xlabel('Julian day');

axB22b=nexttile;
set(gca,'TickDir','out'); 
hold on;

CI=2*std(CUI,0,1,'omitnan');mn=mean(CUI,'omitnan');
% fill([1:365 fliplr(1:365)]',[mn-CI fliplr(mn+CI)]',...
%     [0.9 0.9 0.9],'edgecolor','none','facealpha',0.4);

% plot(1:365,cumsum(BAK.B23,'omitnan'),'color',[0 0 0 0.1],'linewidth',2); 
plot(1:365,CUI,'color',[0.85 0.85 0.85 0.3],'linewidth',1); 
plot(mean(CUI,'omitnan'),'k--','LineWidth',1);
plot(1:365,cumsum(BAK.B22,'omitnan'),'color','k'); 
axis tight; grid on;
% yticks(-400:100:300); 
ylim([-10000 0]);
axis tight
xlim([0 365]);
xticks([0 182 365]);
xlabel('Julian day');
ylabel({'CUI';'(m^3 s^{-1} 100m^{-1})'});
axB22b.YAxis.Exponent = 0;

axB23b=nexttile;
set(gca,'TickDir','out'); 
hold on;
% fill([1:365 fliplr(1:365)]',[mn-CI fliplr(mn+CI)]',...
%     [0.9 0.9 0.9],'edgecolor','none','facealpha',0.4);
plot(1:365,CUI,'color',[0.85 0.85 0.85 0.3],'linewidth',1); 
plot(mean(CUI,'omitnan'),'k--','LineWidth',1);
% plot(1:365,cumsum(BAK.B22,'omitnan'),'color',[0 0 0 0.1],'linewidth',2); 
plot(1:365,cumsum(BAK.B23,'omitnan'),'color','k'); 
axis tight; grid on;
% yticks(-400:100:300); 
ylim([-10000 0]);
axis tight
xlim([0 365]);
xticks([0 182 365]);
%set(gca,'XAxisLocation','top');
% ylabel({'CUI';'(m^3 s^{-1} 100m^{-1})'});

% text(310,155,'2023','fontsize',6,...
%     'HorizontalAlignment','center');
xlabel('Julian day');
axB23b.YAxis.Exponent = 0;

%%%%%%%%%%
% Define labels
labels = {'a)', 'b)', 'c)', 'd)'};

% Define an array to store the axes handles
axHandles = {axB22 axB23 axB22b axB23b};


set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

% Position text labels on each subplot
for idx = 1:length(axHandles)
    axes(axHandles{idx}); % Switch to the current axes

    if idx<3
        lll(1)=line([32 32],[-465 -465+20],'color','k');
        lll(2)=line([60 60],[-465 -465+20],'color','k');
        lll(3)=line([91 91],[-465 -465+20],'color','k');
        lll(4)=line([121 121],[-465 -465+20],'color','k');
        lll(5)=line([152 152],[-465 -465+20],'color','k');
        lll(6)=line([182 182],[-465 -465+20],'color','k');
        lll(7)=line([213 213],[-465 -465+20],'color','k');
        lll(8)=line([244 244],[-465 -465+20],'color','k');
        lll(9)=line([274 274],[-465 -465+20],'color','k');
        lll(10)=line([305 305],[-465 -465+20],'color','k');
        lll(11)=line([335 335],[-465 -465+20],'color','k');

        text([32 60 91 121 152 182 213 244 274 305 335],[-400 -400 -400 -400 -400 -400 -400 -400 -400 -400 -400],...
            {'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'},'fontsize',6,...
            'HorizontalAlignment','center');

        text(0.1, 0.95, labels{idx}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8);

        if idx==1
            text(0.9,0.95,'2022','fontsize',8,'units','normalized',...
            'HorizontalAlignment','center');
        elseif idx==2
            text(0.9,0.95,'2023','fontsize',8,'units','normalized',...
            'HorizontalAlignment','center');
        end


    else
        ll(1)=line([32 32],[-22800 -22700+400],'color','k');
        ll(2)=line([60 60],[-22800 -22700+400],'color','k');
        ll(3)=line([91 91],[-22800 -22700+400],'color','k');
        ll(4)=line([121 121],[-22800 -22700+400],'color','k');
        ll(5)=line([152 152],[-22800 -22700+400],'color','k');
        ll(6)=line([182 182],[-22800 -22700+400],'color','k');
        ll(7)=line([213 213],[-22800 -22700+400],'color','k');
        ll(8)=line([244 244],[-22800 -22700+400],'color','k');
        ll(9)=line([274 274],[-22800 -22700+400],'color','k');
        ll(10)=line([305 305],[-22800 -22700+400],'color','k');
        ll(11)=line([335 335],[-22800 -22700+400],'color','k');
        text([32 60 91 121 152 182 213 244 274 305 335],[-21000 -21000 -21000 -21000 -21000 -21000 -21000 -21000 -21000 -21000 -21000],...
            {'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'},'fontsize',6,...
            'HorizontalAlignment','center');

        text(0.1, 0.95, labels{idx}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8);

        if idx==3
            text(0.9,0.95,'2022','fontsize',8,'units','normalized',...
            'HorizontalAlignment','center');
        elseif idx==4
            text(0.9,0.95,'2023','fontsize',8,'units','normalized',...
            'HorizontalAlignment','center');
        end

        
    end
end



%%%%% Add mooring PP plots
cm1=cmocean('algae');
cm1(1:3,:)=[0.95 0.95 0.95;0.95 0.95 0.95;0.95 0.95 0.95];

msk=M.S2.time>datenum(2023,05,01) & M.S2.time<datenum(2023,08,01) & ...
    ~isnan(M.S2.CTD280m.oxygen) & ~isnan(M.S2.CTD280m.sigma_theta);
Sso=M.S2.CTD280m.oxygen(msk); Sst=M.S2.CTD280m.sigma_theta(msk);
Sao=M.S2.CTD280m.aou(msk);
k = boundary(double(Sso'),double(Sst'));

[N, edgesX, edgesY] = histcounts2(M.S2.CTD280m.oxygen,...
    M.S2.CTD280m.sigma_theta,20:5:250,25:0.05:27.05);
AXh1=nexttile;
set(gca,'TickDir','out'); 
hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p=pcolor(X, Y, N');p.FaceAlpha=0.6;shading flat; colormap(gca,cm1);
xlabel('O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
ylabel('\sigma_\theta (kg m^{-3})');
axis tight;
caxis([0 1300]);
ylim([25 27]);
box on
axis ij
yticks([25:0.5:27])
xlim([20 210]);
xticks(50:50:200);
set(gca,'XTickLabelRotation',0);
text(0.05,0.95,'e)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  
plot(Sso(k), Sst(k), '-', 'LineWidth', 1,'color',rgb_x('violet'));  % Draw boundary
title('Scott2')

[N, edgesX, edgesY] = histcounts2(M.H1.CTD133m.oxygen,...
    M.H1.CTD133m.sigma_theta,20:5:250,25:0.05:27.05);
Hso(isnan(Hso))=[];Hst(isnan(Hst))=[];Ha(isnan(Ha))=[]; 
Ht(isnan(Ht))=[];Hs(isnan(Hs))=[];

% Compute the boundary of the points
k = boundary(double(Hso'),double(Hst'));

% Modify the plotting commands for the oxygen vs sigma theta histogram
AXh2=nexttile;
set(gca, 'TickDir', 'out'); hold on;
X = edgesX(1:end-1) + diff(edgesX)/2;
Y = edgesY(1:end-1) + diff(edgesY)/2;
p = pcolor(X, Y, N'); p.FaceAlpha = 0.6; shading flat; colormap(gca, cm1);
xlabel('O$_2$ ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
% ylabel('\sigma_\theta (kg m^{-3})');
ylim([25 27]); xlim([20 210]);
plot(Hso(k), Hst(k), '-', 'LineWidth', 1,'color',rgb_x('violet'));  % Draw boundary
caxis([0 1300]);
box on;
axis ij;
yticks(25:0.5:27);
xticks(50:50:200);
set(gca, 'XTickLabelRotation', 0);
text(0.05,0.95,'f)','units','normalized');
set(gca, 'Layer', 'top','GridLineStyle', ':');
grid on;  
title('Hak1')


% axes(AXh1)
% c=colorbar('position',[0.425,0.1,0.01,0.075],'fontsize',6);
% c.Ticks=[0 1000];
% title(c,'N','fontsize',6);

axes(AXh2)
c=colorbar('position',[0.85,0.1,0.01,0.075],'fontsize',6);
c.Ticks=[0 1000];
title(c,'N','fontsize',6);

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/BakunUpwell.pdf -dpdf -nofontswap

% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/Bakun.pdf', ...
%                'ContentType', 'vector', 'Resolution', 400);


%% Relation between sigma and o2

msk=~isnan(M.S2.CTD280m.oxygen) & ~isnan(M.S2.CTD280m.sigma_theta);
[rS,p]=corrcoef(double(M.S2.CTD280m.oxygen(msk)),double(M.S2.CTD280m.sigma_theta(msk)));
slopeS=polyfit(double(M.S2.CTD280m.sigma_theta(msk))',double(M.S2.CTD280m.oxygen(msk))',1);


msk=~isnan(M.H1.CTD133m.oxygen) & ~isnan(M.H1.CTD133m.sigma_theta);
[rH,p]=corrcoef(double(M.H1.CTD133m.oxygen(msk)),double(M.H1.CTD133m.sigma_theta(msk)));
slopeH=polyfit(double(M.H1.CTD133m.sigma_theta(msk))',double(M.H1.CTD133m.oxygen(msk))',1);

%% what about STI and TUMI? ALGORITHM

load isoDO_DFOandHakai.mat
sig_thresh = [26 26.5 26.75];
isoDO.oxygen(isoDO.oxygen < 15) = isoDO.oxygen(isoDO.oxygen < 15) * 43.7;
yrs = 2003:2023;
num_years = length(yrs);
num_sig_thresh = length(sig_thresh);
isoDO.oxygen(isoDO.oxygen>600)=NaN;

% Initialize matrices for annual means and 95% CI for depth, temperature, and salinity
isoDO.AmnDO = NaN(num_sig_thresh, num_years);

for i = 2003:2023
    year_idx = i - 2002;
    msk = isoDO.time > datenum(i, 1, 1) & isoDO.time < datenum(i + 1, 1, 1) & isoDO.latitude < 52 & ...
        isoDO.OSdist>80;

    for j = 1:num_sig_thresh
        % Dissolved Oxygen
        data_DO = isoDO.oxygen(j, msk);
        data_DO = data_DO(~isnan(data_DO)); % Remove NaNs
        if ~isempty(data_DO)
            isoDO.AmnDO(j, year_idx) = mean(data_DO);
        end 
    end
end

% Adjust the smoothing to a 15-day running average
smthCUI = smoothdata(CUI(37:end,:), 2, 'movmean', 15);

STI = NaN(21,1); % Initialize the Spring Transition Index array
TUMI = NaN(21,1); % Initialize the Spring Transition Index array
FTI = NaN(21,1); % Initialize the Spring Transition Index array

lnes=linspecer(21);

figure; hold on;
for i = 1:21
    gradientCUI = diff(smthCUI(i,1:end)); % Calculate the gradient
    for j = 1:length(gradientCUI) - 15 % Adjust loop to check for the condition
        if gradientCUI(j) > 0 && all(gradientCUI(j+1:j+15) > 0) % Check if the gradient is positive for at least 14 days
            STI(i) = j+0; % Store the day of transition
            break; % Stop the loop once the condition is met
        end
    end

    for j = STI(i):length(gradientCUI) - 15 % Adjust loop to check for the condition
        if gradientCUI(j) < 0 && all(gradientCUI(j+1:j+15) < 0) % Check if the gradient is positive for at least 14 days
            FTI(i) = j; % Store the day of transition
            break; % Stop the loop once the condition is met
        end
    end
    try
        
        [~,idx]=max(CUI(36+i,STI(i):end));

        TUMI(i,1)=sum(smthCUI(STI(i):FTI(i)),'omitnan');
        plot(1:365, smthCUI(i,:), 'LineWidth', 1.5,'Color',lnes(i,:)'); % Plot the smoothed CUI
        if ~isnan(STI(i))
            scatter(STI(i), smthCUI(i,STI(i)), 40, 'filled', 'd','MarkerFaceColor',lnes(i,:)); % Mark the transition
            scatter(FTI(i), smthCUI(i,FTI(i)), 40, 'filled', 's','MarkerFaceColor',lnes(i,:)); % Mark the transition

        end
    end
end
clc
fprintf('STI=%3.0f+-%2.0f days\n',nanmean(STI),nanstd(STI));
fprintf('FTI=%3.0f+-%2.0f days\n',nanmean(FTI),nanstd(FTI));

xlabel('Day of the Year');
ylabel('Smoothed CUI');
title('Spring Transition of Upwelling');
legend('Smoothed CUI', 'Transition Point');
grid on;

figure
lnes=linspecer(3);
for i=1:3
    subplot(3,3,i); hold on;
    scatter(STI,isoDO.AmnDO(i,:),'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:)) & ~isnan(STI)';
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),STI(msk));
    text(0.3,0.1,[num2str(r(2)),'/',num2str(pval(2))],'units','normalized')
    title('STI vs. ~O')
    ylabel('O2');
xlabel('STI');
end

tmp=STI;tmp(STI<45)=NaN;
lnes=linspecer(3);
for i=1:3
    subplot(3,3,i+3); hold on;
    scatter(tmp,isoDO.AmnDO(i,:),'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:)) & ~isnan(tmp)';
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),tmp(msk));
    text(0.3,0.1,[num2str(r(2)),'/',num2str(pval(2))],'units','normalized')
    title('STI (early removed) vs. ~O')
    ylabel('O2');
xlabel('STI');
end
fprintf('STI (early removed)=%3.0f+-%2.0f days\n',nanmean(tmp),nanstd(tmp));

lnes=linspecer(3);
for i=1:3
    subplot(3,3,i+6); hold on;
    scatter(TUMI,isoDO.AmnDO(i,:),'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:)) & ~isnan(TUMI)';
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),TUMI(msk));
    text(0.3,0.1,[num2str(r(2)),'/',num2str(pval(2))],'units','normalized')
    title('TUMI vs. ~O')

end

%% Figure;
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 15 10],'color','w');
lnes=linspecer(3);
for i=1:3
    subplot(1,3,i); hold on;
    scatter(STI,isoDO.AmnDO(i,:),'filled','markerfacecolor',lnes(i,:));
    msk=~isnan(isoDO.AmnDO(i,:)) & ~isnan(STI)';
    [r,pval]=corrcoef(isoDO.AmnDO(i,msk),STI(msk));
    text(0.45,0.9,{['r=',num2str(r(2),'%2.2f')];...
        ['p-value=',num2str(pval(2),'%2.2f')]},'units','normalized',...
        'horizontalalignment','right');
    title([num2str(sig_thresh(i),'%2.2f') ' \sigma_\theta']);
    ylabel('$\overline{\rm{O}_2}$ ($\mu$mol~kg$^{-1}$)','interpreter','latex');
    grid on;
    xlabel('STI (day of year)');
end

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);
%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/STIvsO2.jpg -m4 -nofontswap

%% Correlate NPGO
load aNPGO

f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 8 8],'color','w');
lnes=linspecer(3);
scatter(STI,aNPGO,'filled','markerfacecolor','k');
msk=~isnan(aNPGO) & ~isnan(STI);
[r,pval]=corrcoef(aNPGO(msk),STI(msk));
text(0.45,0.9,{['r=',num2str(r(2),'%2.2f')];...
    ['p-value=',num2str(pval(2),'%2.2f')]},'units','normalized',...
    'horizontalalignment','right');
ylabel('NPGO Index');
grid on;
xlabel('STI (day of year)');

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/STIvsNPGO.jpg -m4 -nofontswap


%% sensitivity test
mnSTI=[];
for k=2:30
    STI = NaN(21,1); % Initialize the Spring Transition Index array
    TUMI = NaN(21,1); % Initialize the Spring Transition Index array

    for i = 1:21
        gradientCUI = diff(smthCUI(i,1:end)); % Calculate the gradient
        for j = 1:length(gradientCUI) - k % Adjust loop to check for the condition
            if gradientCUI(j) > 0 && all(gradientCUI(j+1:j+k) > 0) % Check if the gradient is positive for at least 14 days
                STI(i) = j; % Store the day of transition
                break; % Stop the loop once the condition is met
            end
        end
        try
            STI(i)=STI(i);
            [~,idx]=max(CUI(36+i,STI(i):end));

            TUMI(i,1)=sum(CUI(36+i,STI(i):idx),'omitnan');
        end
    end
    mnSTI(k-1)=nanmean(STI);
end
figure; plot(2:30,mnSTI)
