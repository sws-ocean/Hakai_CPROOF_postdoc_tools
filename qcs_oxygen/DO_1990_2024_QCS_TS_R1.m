%% Load in all QCS oxygen

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear
load('Dcorrected_allProfs_20241011.mat');

%% Load bathymetry

lon_lim=[-132 -127];
lat_lim=[50 53];
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1])';

[l,ll]=meshgrid(lon(ilon),lat(ilat));

% Offshore-axis
osl=linspace(-129.3360,-127.5683,1000);
osll=linspace(52.7003,51.1757,1000);

%% Run through file list and check for oxygen names

folder = '/Users/samuelstevens/Dropbox/Hakai/ctd/data/DO_1990-2024_QCS/outdir/';
files = dir(fullfile(folder, '*.nc'));  % List all .nc files in the directory

oxygen_vars = {};  % Cell array to store names of oxygen-related variables

for i = 1:length(files)
    filename = fullfile(folder, files(i).name);
    info = ncinfo(filename);

    for j = 1:length(info.Variables)
        varName = info.Variables(j).Name;
        if contains(varName, 'O2', 'IgnoreCase', true) || ...
                contains(varName, 'ox', 'IgnoreCase', true) || ...
                contains(varName, 'DO', 'IgnoreCase', true)
            oxygen_vars{end+1} = varName;  % Add to list if it's an oxygen-related variable
        end
    end
    % ncdisp(filename)
    % pause
    % clc
end

oxygen_vars = unique(oxygen_vars);  % Remove duplicates
disp(oxygen_vars); % all o2 vars contain OX, DOXMZZ01 is umol/kg


%% Load data
folder = '/Users/samuelstevens/Dropbox/Hakai/ctd/data/DO_1990-2024_QCS/outdir/';
files = dir(fullfile(folder, '*.nc'));  % List all .nc files

% Initialize a structure to store the data
data = struct();
sig_thresh=[26 26.5 26.75];

isoDO=struct();

isoDO.oxygen=NaN(length(sig_thresh),length(files));
isoDO.rho=NaN(length(sig_thresh),length(files));
isoDO.z=NaN(length(sig_thresh),length(files));
isoDO.t=NaN(length(sig_thresh),length(files));
isoDO.SA=NaN(length(sig_thresh),length(files));
isoDO.AOU=NaN(length(sig_thresh),length(files));
isoDO.unit=cell(1,length(files));
isoDO.longitude=NaN(1,length(files));
isoDO.latitude=NaN(1,length(files));
isoDO.time=NaN(1,length(files));
isoDO.h=NaN(1,length(files));
isoDO.Poxygen=NaN(350,length(files));
isoDO.Pz=repmat([1:350]',1,length(files));

bottomDO=struct();

h = waitbar(0, 'Processing files...');  % Initialize the waitbar
count=0;
count2=0;
for i = 1:length(files)
    filename = fullfile(folder, files(i).name);
    info = ncinfo(filename);
    waitbar(i / length(files), h, sprintf('Processing %d of %d files...', i, length(files)));


    % Create a valid MATLAB field name: replace invalid characters and ensure it starts with a letter
    validFieldName = regexprep(files(i).name, '[^a-zA-Z0-9_]', '_');
    if ~isletter(validFieldName(1))
        validFieldName = ['Data_' validFieldName];  % Prefix with 'Data_' if not starting with a letter
    end

    % Temporary structure for each file's data
    file_data = struct();

    for j = 1:length(info.Variables)
        varName = info.Variables(j).Name;

        % Read the variable data
        file_data.(varName) = ncread(filename, varName);
    end

    try
        % Oxygen fieldnames
        flds=fieldnames(file_data);
        idxO=find(contains(flds,'ox','IgnoreCase', true));

        if sum(contains({flds{idxO}},'oxm','IgnoreCase', true))
            idxO=find(contains(flds,'oxm','IgnoreCase', true));
            unit='umol';
        elseif isscalar(idxO) && contains({flds{idxO}},'oxy','IgnoreCase', true)
            unit='ml/l';
        else
            idxO=idxO(1);
            unit='ml/l';
        end

        % find bottom depth
        hh=interp2(l,ll,Z,file_data.longitude,file_data.latitude,'nearest');
    
        % find bottom depth
        hh=interp2(l,ll,Z,file_data.longitude,file_data.latitude,'nearest');

        if abs(hh)-max(abs(file_data.depth))<20 &&...
                contains(file_data.instrument_type,'Bottle','IgnoreCase', true)
            count=count+1;
            bottomDO.h(count)=hh;
            bottomDO.unit{count}=unit;
            [bottomDO.z(count),idxZ]=max(file_data.depth);
            bottomDO.oxygen(count)=file_data.(flds{idxO(1)})(idxZ);
            bottomDO.longitude(count)=file_data.longitude;
            bottomDO.latitude(count)=file_data.latitude;
            bottomDO.time(count)=file_data.time;

            try
                SA=gsw_SA_from_SP(file_data.sea_water_practical_salinity,...
                    file_data.sea_water_pressure,file_data.longitude,...
                    file_data.latitude);
                CT=gsw_CT_from_t(SA,file_data.sea_water_temperature,...
                    file_data.sea_water_pressure);
                sigth=gsw_sigma0(SA,CT);
                bottomDO.s_t(count)=max(sigth);
            end

        elseif abs(hh)-max(abs(file_data.depth))<20 &&...
                contains(file_data.instrument_type,'CTD','IgnoreCase', true)
            count=count+1;
            bottomDO.h(count)=hh;
            bottomDO.unit{count}=unit;
            bottomDO.z(count)=max(file_data.depth);
            msk=(abs(hh)-abs(file_data.depth))<20;
            bottomDO.oxygen(count)=min(file_data.(flds{idxO(1)})(msk));
            bottomDO.longitude(count)=file_data.longitude;
            bottomDO.latitude(count)=file_data.latitude;
            bottomDO.time(count)=file_data.time;
            
            try
                SA=gsw_SA_from_SP(file_data.sea_water_practical_salinity,...
                    file_data.sea_water_pressure,file_data.longitude,...
                    file_data.latitude);
                CT=gsw_CT_from_t(SA,file_data.sea_water_temperature,...
                    file_data.sea_water_pressure);
                sigth=gsw_sigma0(SA,CT);
                bottomDO.s_t(count)=max(sigth);
            end
        end


        try
            SA=gsw_SA_from_SP(file_data.sea_water_practical_salinity,...
                file_data.sea_water_pressure,file_data.longitude,...
                file_data.latitude);
            CT=gsw_CT_from_t(SA,file_data.sea_water_temperature,...
                file_data.sea_water_pressure);
            sigth=gsw_sigma0(SA,CT);
            msk=~isnan(sigth);

            if max(sigth>=sig_thresh(1)) && hh>-350

                rho=gsw_rho(SA,CT,file_data.sea_water_pressure);
                PT=gsw_pt0_from_t(SA,file_data.sea_water_temperature,...
                    file_data.sea_water_pressure);
                Caou=aou(file_data.sea_water_practical_salinity,...
                    PT,file_data.(flds{idxO(1)}));

                count2=count2+1;
                isoDO.oxygen(:,count2)=interp1(sigth(msk),file_data.(flds{idxO(1)})(msk),...
                    sig_thresh);
                isoDO.z(:,count2)=interp1(sigth(msk),file_data.depth(msk),sig_thresh);
                isoDO.t(:,count2)=interp1(sigth(msk),file_data.sea_water_temperature(msk),...
                    sig_thresh);
                isoDO.AOU(:,count2)=interp1(sigth(msk),Caou(msk),...
                    sig_thresh);
                isoDO.SA(:,count2)=interp1(sigth(msk),SA(msk),sig_thresh);
                isoDO.rho(:,count2)=interp1(sigth(msk),rho(msk),sig_thresh);
                isoDO.unit{count2}=unit;
                isoDO.longitude(count2)=file_data.longitude;
                isoDO.latitude(count2)=file_data.latitude;
                isoDO.time(count2)=file_data.time;
                isoDO.h(count2)=hh;
                try
                isoDO.Poxygen(:,count2)=interp1(file_data.depth(msk)+sort((rand(size(file_data.depth(msk)))/1000)),...
                    file_data.(flds{idxO(1)})(msk),[1:350]');
                catch
                    keyboard
                end
            end
        catch
            warning('Could not calculate isopycnal properties for %s!',...
                filename)

        end
    end

end

close(h);  % Close the waitbar

mt=datetime(bottomDO.time,'ConvertFrom', 'posix', 'TimeZone', 'UTC');
bottomDO.time=datenum(mt);
msk=strcmp(bottomDO.unit,'ml/l');
bottomDO.oxygen(msk)=bottomDO.oxygen(msk).*43.7;
bottomDO.s_t(bottomDO.s_t==0)=NaN;

msk=bottomDO.oxygen<0;
bottomDO.h(msk)=[];
bottomDO.z(msk)=[];
bottomDO.oxygen(msk)=[];
bottomDO.longitude(msk)=[];
bottomDO.latitude(msk)=[];
bottomDO.time(msk)=[];
bottomDO.s_t(msk)=[];

mt=datetime(isoDO.time,'ConvertFrom', 'posix', 'TimeZone', 'UTC');
isoDO.time=datenum(mt);

msk=strcmp(isoDO.unit,'ml/l');
isoDO.oxygen(msk) = isoDO.oxygen(msk) .* (44.661 ./ isoDO.rho(msk));
isoDO.Poxygen(:,msk) = isoDO.Poxygen(:,msk) .* (44.661 ./ isoDO.rho(msk));

msk=isoDO.oxygen<0;
isoDO.oxygen(msk)=NaN;
isoDO.t(msk)=NaN;
isoDO.SA(msk)=NaN;

msk=isnan(isoDO.time) | isoDO.time<datenum(2003,01,01);
isoDO.h(msk)=[];
isoDO.z(:,msk)=[];
isoDO.oxygen(:,msk)=[];
isoDO.longitude(msk)=[];
isoDO.latitude(msk)=[];
isoDO.time(msk)=[];
isoDO.t(:,msk)=[];
isoDO.SA(:,msk)=[];
isoDO.rho(:,msk)=[];
isoDO.AOU(:,msk)=[];
isoDO.Poxygen(:,msk)=[];
isoDO.Pz(:,msk)=[];

isoDO.oxygen(isoDO.oxygen>1000)=isoDO.oxygen(isoDO.oxygen>1000)/43.7;
isoDO.oxygen(isoDO.oxygen<15)=isoDO.oxygen(isoDO.oxygen<15)*43.7;
isoDO.oxygen(isoDO.oxygen>250)=NaN;
isoDO.oxygen(2,isoDO.oxygen(2,:)>160)=NaN;

l=[ -127.7593 -129.343032262174	-129.778110234244	-130.352413157377	-130.665669297267	-130.178381968549	-129.551869688768	-129.169001073346	-128.420666961386	-128.072604583730	-127.620123492777	-128.159620178144	-128.768729339042];
ll=[50.9022 52.6069254185693	52.3725266362253	52.0742009132420	51.7332572298326	51.1792237442922	50.6997716894977	50.6571537290715	50.9341704718417	51.1792237442922	51.6160578386606	52.1914003044140	52.6069254185693];
AS=alphaShape(l',ll',1,'HoleThreshold',15);

msk=~inShape(AS,isoDO.longitude,isoDO.latitude);
isoDO.h(msk)=[];
isoDO.z(:,msk)=[];
isoDO.oxygen(:,msk)=[];
isoDO.longitude(msk)=[];
isoDO.latitude(msk)=[];
isoDO.time(msk)=[];
isoDO.t(:,msk)=[];
isoDO.SA(:,msk)=[];
isoDO.rho(:,msk)=[];
isoDO.AOU(:,msk)=[];
isoDO.Poxygen(:,msk)=[];
isoDO.Pz(:,msk)=[];


bottomDO.oxygen(bottomDO.oxygen>1000)=bottomDO.oxygen(bottomDO.oxygen>1000)/43.7;
bottomDO.oxygen(bottomDO.oxygen<15)=bottomDO.oxygen(bottomDO.oxygen<15)*43.7;
msk=~inShape(AS,bottomDO.longitude,bottomDO.latitude);
bottomDO.h(msk)=[];
bottomDO.z(msk)=[];
bottomDO.oxygen(msk)=[];
bottomDO.longitude(msk)=[];
bottomDO.latitude(msk)=[];
bottomDO.time(msk)=[];
bottomDO.s_t(msk)=[];

% add offshore distance variable
for i=1:length(isoDO.AOU)
    d=gsw_distance([ones(size(osl)).*isoDO.longitude(i); osl]',...
        [ones(size(osl)).*isoDO.latitude(i); osll]');
    isoDO.OSdist(i)=min(d);
end
figure;scatter(isoDO.OSdist,isoDO.AOU(2,:));
figure;scatter(isoDO.time,isoDO.oxygen(2,:));
figure; plot(isoDO.Poxygen,isoDO.Pz); axis ij; grid on;
% save('/Users/samuelstevens/Dropbox/Hakai/ctd/data/DO_1990-2024_QCS/DO_1990-2024_QCS_v3.mat','isoDO','bottomDO');

%% Add Hakai data
load DO_1990-2024_QCS_v3.mat
load Hakai_allQCSPLUS.mat
sig_thresh = [26 26.5 26.75];

lon_lim=[-132 -127];
lat_lim=[50 53];
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1])';

[l,ll]=meshgrid(lon(ilon),lat(ilat));

% Offshore-axis
osl=linspace(-129.3360,-127.5683,1000);
osll=linspace(52.7003,51.1757,1000);

H = HakaiWaterPropertiesInstrumentP;
clearvars HakaiWaterPropertiesInstrumentP
H.hakai_id = string(H.hakai_id);
flds = fieldnames(H);

H.time = datenum(1970,01,01,00,00,00) + (H.time / 86400);
[H.allT, idx] = unique(H.time);
H.l = H.longitude(idx);
H.ll = H.latitude(idx);

isoH.oxygen=NaN(length(sig_thresh),length(H.allT));
isoH.z=NaN(length(sig_thresh),length(H.allT));
isoH.t=NaN(length(sig_thresh),length(H.allT));
isoH.longitude=NaN(1,length(H.allT));
isoH.latitude=NaN(1,length(H.allT));
isoH.time=NaN(1,length(H.allT));
isoH.AOU=NaN(length(sig_thresh),length(H.allT));

bottomH.oxygen=NaN(1,length(H.allT));
bottomH.z=NaN(1,length(H.allT));
bottomH.h=NaN(1,length(H.allT));
bottomH.t=NaN(1,length(H.allT));
bottomH.longitude=NaN(1,length(H.allT));
bottomH.latitude=NaN(1,length(H.allT));
bottomH.time=NaN(1,length(H.allT));
bottomH.s_t=NaN(1,length(H.allT));

H.pressure=double(H.pressure);

H.SA=gsw_SA_from_SP(H.salinity,H.pressure,H.longitude,H.latitude);
H.CT=gsw_CT_from_t(H.SA,H.temperature,H.pressure);
H.rho=gsw_rho(H.SA,H.CT,H.pressure);
H.DO_umol_kg = (H.dissolved_oxygen_ml_l * 44660) ./ H.rho;
H.s_t=gsw_sigma0(H.SA,H.CT);
H.pt=gsw_pt_from_t(H.SA,H.temperature,H.pressure);
H.AOU=aou(H.salinity,H.pt,H.dissolved_oxygen_ml_l);

count=0;
count2=0;

h = waitbar(0, 'Processing files...');  % Initialize the waitbar
for i = 1:length(H.l)
    waitbar(i / length(H.l), h, sprintf('Processing %d of %d files...', i, length(H.l)));

    cmsk = H.time == H.allT(i);
    
    tmpDO=H.DO_umol_kg(cmsk);
    tmpST=H.s_t(cmsk);
    tmpT=H.temperature(cmsk);
    tmpZ=H.depth(cmsk);
    tmpAOU=H.AOU(cmsk);

    if sum(~isnan(tmpDO))>2 & length(tmpDO)>1
        count=count+1;
        isoH.oxygen(:,count)=interp1(tmpST(~isnan(tmpDO))+(rand(size(tmpST(~isnan(tmpDO))))./1000),...
            tmpDO(~isnan(tmpDO)),sig_thresh);
        isoH.t(:,count)=interp1(tmpST(~isnan(tmpDO))+(rand(size(tmpST(~isnan(tmpDO))))./1000),...
            tmpT(~isnan(tmpDO)),sig_thresh);
        isoH.z(:,count)=interp1(tmpST(~isnan(tmpDO))+(rand(size(tmpST(~isnan(tmpDO))))./1000),...
            tmpZ(~isnan(tmpDO)),sig_thresh);
        isoH.longitude(count)=mean(H.longitude(cmsk));
        isoH.latitude(count)=mean(H.latitude(cmsk));
        isoH.time(count)=mean(H.allT(i));
        isoH.AOU(:,count)=interp1(tmpST(~isnan(tmpDO))+(rand(size(tmpST(~isnan(tmpDO))))./1000),...
            tmpAOU(~isnan(tmpDO)),sig_thresh);
    end

    if sum(~isnan(tmpDO))>2 && length(tmpDO)>1 && max(tmpZ)>75 &&...
            min(tmpZ)<15

        hh=interp2(l,ll,Z,mean(H.longitude(cmsk)),mean(H.latitude(cmsk)),'nearest');

        if abs(hh)-max(abs(tmpZ))<20
            count2=count2+1;
            bottomH.oxygen(:,count2)=min(tmpDO(end-20:end));
            bottomH.z(:,count2)=max(tmpZ);
            bottomH.longitude(count2)=mean(H.longitude(cmsk));
            bottomH.latitude(count2)=mean(H.latitude(cmsk));
            bottomH.time(count2)=mean(H.allT(i));
            bottomH.h(count2)=hh;
            bottomH.s_t(count2)=max(tmpST);
        end
    end
end
close(h);  % Close the waitbar

isoH.t(:,sum(~isnan(isoH.oxygen))==0)=[];
isoH.z(:,sum(~isnan(isoH.oxygen))==0)=[];
isoH.time(sum(~isnan(isoH.oxygen))==0)=[];
isoH.AOU(:,sum(~isnan(isoH.oxygen))==0)=[];
isoH.longitude(sum(~isnan(isoH.oxygen))==0)=[];
isoH.latitude(sum(~isnan(isoH.oxygen))==0)=[];
isoH.oxygen(:,sum(~isnan(isoH.oxygen))==0)=[];


% add offshore distance variable
for i=1:length(isoH.oxygen)
    d=gsw_distance([ones(size(osl)).*isoH.longitude(i); osl]',...
        [ones(size(osl)).*isoH.latitude(i); osll]');
    isoH.OSdist(i)=min(d);
end


% add H to isoDO dataset
isoDO.oxygen=[isoDO.oxygen isoH.oxygen];
isoDO.t=[isoDO.t isoH.t];
isoDO.z=[isoDO.z isoH.z];
isoDO.time=[isoDO.time isoH.time];
isoDO.longitude=[isoDO.longitude isoH.longitude];
isoDO.latitude=[isoDO.latitude isoH.latitude];
isoDO.OSdist=[isoDO.OSdist isoH.OSdist];
isoDO.AOU=[isoDO.AOU isoH.AOU];

bottomH.t(:,isnan(bottomH.oxygen))=[];
bottomH.z(:,isnan(bottomH.oxygen))=[];
bottomH.h(:,isnan(bottomH.oxygen))=[];
bottomH.time(:,isnan(bottomH.oxygen))=[];
bottomH.longitude(:,isnan(bottomH.oxygen))=[];
bottomH.latitude(:,isnan(bottomH.oxygen))=[];
bottomH.s_t(:,isnan(bottomH.oxygen))=[];
bottomH.oxygen(:,isnan(bottomH.oxygen))=[];

% add H to bottomDO dataset
bottomDO.oxygen=[bottomDO.oxygen bottomH.oxygen];
bottomDO.z=[bottomDO.z bottomH.z];
bottomDO.h=[bottomDO.h bottomH.h];
bottomDO.time=[bottomDO.time bottomH.time];
bottomDO.longitude=[bottomDO.longitude bottomH.longitude];
bottomDO.latitude=[bottomDO.latitude bottomH.latitude];
bottomDO.s_t=[bottomDO.s_t bottomH.s_t];

% Remove winter values
mnth=month(isoDO.time);
isoDO.oxygen(:,mnth<5 | mnth>10)=NaN;
isoDO.t(:,mnth<5 | mnth>10)=NaN;
isoDO.z(:,mnth<5 | mnth>10)=NaN;

% save('/Users/samuelstevens/Dropbox/Hakai/gliders/data/isoDO_DFOandHakai.mat','isoDO','bottomDO');

%% Isopycnal Time series analysis
load isoDO_DFOandHakai.mat
sig_thresh = [26 26.5 26.75];
isoDO.oxygen(isoDO.oxygen < 15) = isoDO.oxygen(isoDO.oxygen < 15) * 43.7;
isoDO.oxygen(isoDO.oxygen > 700) = NaN;

yrs = 2003:2023;
num_years = length(yrs);
num_sig_thresh = length(sig_thresh);

% clc

% Initialize matrices for annual means and 95% CI for depth, temperature, and salinity
isoDO.AmnDO = NaN(num_sig_thresh, num_years);
isoDO.ACI_DO = NaN(num_sig_thresh, 2, num_years); % Confidence intervals for DO
isoDO.AmnDepth = NaN(num_sig_thresh, num_years);
isoDO.ACI_Depth = NaN(num_sig_thresh, 2, num_years); % Confidence intervals for Depth
isoDO.AmnTemp = NaN(num_sig_thresh, num_years);
isoDO.ACI_Temp = NaN(num_sig_thresh, 2, num_years); % Confidence intervals for Temp
isoDO.AmnSal = NaN(num_sig_thresh, num_years);
isoDO.ACI_Sal = NaN(num_sig_thresh, 2, num_years); % Confidence intervals for Salinity
isoDO.AmnAOU = NaN(num_sig_thresh, num_years);
isoDO.ACI_AOU = NaN(num_sig_thresh, 2, num_years); % Confidence intervals for AOU

nBoot = 1000; % Number of bootstrap samples
for i = 2003:2023
    year_idx = i - 2002;
    msk = isoDO.time > datenum(i, 1, 1) & isoDO.time < datenum(i + 1, 1, 1);

    for j = 1:num_sig_thresh
        % Dissolved Oxygen
        data_DO = isoDO.oxygen(j, msk);
        data_DO = data_DO(~isnan(data_DO)); % Remove NaNs

        if ~isempty(data_DO)
            isoDO.AmnDO(j, year_idx) = mean(data_DO);
            if length(data_DO) > 1
                bootstat_DO = bootstrp(nBoot, @mean, data_DO);
                isoDO.ACI_DO(j, :, year_idx) = prctile(bootstat_DO, [2.5 97.5]);
                % isoDO.AseDO(j, year_idx)=diff(prctile(bootstat_DO, [2.5 97.5]))/2;
                isoDO.AseDO(j, year_idx)=nanstd(data_DO);
            else
                isoDO.ACI_DO(j, :, year_idx) = [NaN NaN]; % Cannot compute CI with one data point
                isoDO.AseDO(j, year_idx)=0;
            end
        end

        % Depth
        data_Depth = isoDO.z(j, msk);
        data_Depth = data_Depth(~isnan(data_Depth));
        if ~isempty(data_Depth)
            isoDO.AmnDepth(j, year_idx) = mean(data_Depth);
            if length(data_Depth) > 1
                bootstat_Depth = bootstrp(nBoot, @mean, data_Depth);
                % isoDO.ACI_Depth(j, :, year_idx) = prctile(bootstat_Depth, [2.5 97.5]);
                isoDO.ACI_Depth(j, :, year_idx) = nanstd(bootstat_Depth);
                
            else
                isoDO.ACI_Depth(j, :, year_idx) = [NaN NaN];
            end
        end

        % Temperature
        data_Temp = isoDO.t(j, msk);
        data_Temp = data_Temp(~isnan(data_Temp));
        if ~isempty(data_Temp)
            isoDO.AmnTemp(j, year_idx) = mean(data_Temp);
            if length(data_Temp) > 1
                bootstat_Temp = bootstrp(nBoot, @mean, data_Temp);
                % isoDO.ACI_Temp(j, :, year_idx) = prctile(bootstat_Temp, [2.5 97.5]);
                isoDO.ACI_Temp(j, :, year_idx) = nanstd(bootstat_Temp);

            else
                isoDO.ACI_Temp(j, :, year_idx) = [NaN NaN];
            end
        end

        % % Salinity
        % data_Sal = isoDO.SA(j, msk);
        % data_Sal = data_Sal(~isnan(data_Sal));
        % if ~isempty(data_Sal)
        %     isoDO.AmnSal(j, year_idx) = mean(data_Sal);
        %     if length(data_Sal) > 1
        %         bootstat_Sal = bootstrp(nBoot, @mean, data_Sal);
        %         isoDO.ACI_Sal(j, :, year_idx) = prctile(bootstat_Sal, [2.5 97.5]);
        %     else
        %         isoDO.ACI_Sal(j, :, year_idx) = [NaN NaN];
        %     end
        % end

        % % Apparent Oxygen Utilization (AOU)
        % data_AOU = isoDO.AOU(j, msk);
        % data_AOU = data_AOU(~isnan(data_AOU));
        % if ~isempty(data_AOU)
        %     isoDO.AmnAOU(j, year_idx) = mean(data_AOU);
        %     if length(data_AOU) > 1
        %         bootstat_AOU = bootstrp(nBoot, @mean, data_AOU);
        %         isoDO.ACI_AOU(j, :, year_idx) = prctile(bootstat_AOU, [2.5 97.5]);
        %     else
        %         isoDO.ACI_AOU(j, :, year_idx) = [NaN NaN];
        %     end
        % end
    end
end

% Initialize matrices for linear trend calculations and confidence intervals
isoDO.DOtrend = NaN(1, num_sig_thresh);
isoDO.DepthTrend = NaN(1, num_sig_thresh);
isoDO.TempTrend = NaN(1, num_sig_thresh);
isoDO.SalTrend = NaN(1, num_sig_thresh);
isoDO.DOtrendCI = NaN(2, num_sig_thresh); % 95% CI for DO trend
isoDO.DepthTrendCI = NaN(2, num_sig_thresh); % 95% CI for Depth trend
isoDO.TempTrendCI = NaN(2, num_sig_thresh); % 95% CI for Temp trend
isoDO.SalTrendCI = NaN(2, num_sig_thresh); % 95% CI for Sal trend
isoDO.fittedy = NaN(length(yrs), num_sig_thresh);

predictYrs=2024:2200;
isoDO.predictY = NaN(length(predictYrs), num_sig_thresh);
isoDO.predictY_low = NaN(length(predictYrs), num_sig_thresh);
isoDO.predictY_high = NaN(length(predictYrs), num_sig_thresh);

% Number of bootstrap samples
nBoot = 1000;

% Calculate linear trends using first-order polynomial (linear regression) and bootstrap for confidence intervals
for i = 1:num_sig_thresh
    % Dissolved Oxygen trend
    x = yrs(~isnan(isoDO.AmnDO(i, :)));
    y = isoDO.AmnDO(i, ~isnan(isoDO.AmnDO(i, :)));
    p = tsreg(x,y);
    isoDO.DOtrend(1, i) = p(1); % Slope of the linear fit    
    bootFun = @(bootr) polyfit(x(bootr), y(bootr), 1);
    bootTrends = bootstrp(nBoot, bootFun, 1:length(x));
    % isoDO.DOtrendCI(:, i) = prctile(bootTrends(:, 1), [2.5 97.5]);
    isoDO.DOtrendCI(:, i) = [isoDO.DOtrend(1, i)-nanstd(bootTrends(:, 1))*2 ...
        isoDO.DOtrend(1, i)+nanstd(bootTrends(:, 1))*2];

    isoDO.fittedy(:,i)=polyval(mean(bootTrends),yrs,1);
    isoDO.predictY(:,i)=polyval(mean(bootTrends),predictYrs,1);
    isoDO.predictY_low(:, i) = isoDO.fittedy(end,i)+(isoDO.DOtrendCI(2, i).*...
        (predictYrs-predictYrs(1)));
    isoDO.predictY_high(:, i) = isoDO.fittedy(end,i)+(isoDO.DOtrendCI(1, i).*...
        (predictYrs-predictYrs(1)));

    [~,isoDO.DOp(1, i)]=mann_kendall(y);

    % Depth trend
    x = yrs(~isnan(isoDO.AmnDepth(i, :)));
    y = isoDO.AmnDepth(i, ~isnan(isoDO.AmnDepth(i, :)));
    p = tsreg(x,y);
    isoDO.DepthTrend(1, i) = p(1); % Slope of the linear fit
    bootFun = @(bootr) polyfit(x(bootr), y(bootr), 1);
    bootTrends = bootstrp(nBoot, bootFun, 1:length(x));
    % isoDO.DepthTrendCI(:, i) = prctile(bootTrends(:, 1), [2.5 97.5]);
    isoDO.DepthTrendCI(:, i) =[isoDO.DepthTrend(1, i)-nanstd(bootTrends(:, 1))*2 ...
        isoDO.DepthTrend(1, i)+nanstd(bootTrends(:, 1))*2];
    [~,isoDO.depthp(1, i)]=mann_kendall(y);
    
    % Temperature trend
    x = yrs(~isnan(isoDO.AmnTemp(i, :)));
    y = isoDO.AmnTemp(i, ~isnan(isoDO.AmnTemp(i, :)));
    p = tsreg(x,y);
    isoDO.TempTrend(1, i) = p(1); % Slope of the linear fit
    bootFun = @(bootr) polyfit(x(bootr), y(bootr), 1);
    bootTrends = bootstrp(nBoot, bootFun, 1:length(x));
    isoDO.TempTrendCI(:, i) =[isoDO.TempTrend(1, i)-nanstd(bootTrends(:, 1))*2 ...
        isoDO.TempTrend(1, i)+nanstd(bootTrends(:, 1))*2];
    % isoDO.TempTrendCI(:, i) = prctile(bootTrends(:, 1), [2.5 97.5]);
    [~,isoDO.tempp(1, i)]=mann_kendall(y);

    
    % % Salinity trend
    % x = yrs(~isnan(isoDO.AmnSal(i, :)));
    % y = isoDO.AmnSal(i, ~isnan(isoDO.AmnSal(i, :)));
    % p = tsreg(x,y);
    % isoDO.SalTrend(1, i) = p(1); % Slope of the linear fit
    % bootFun = @(bootr) polyfit(x(bootr), y(bootr), 1);
    % bootTrends = bootstrp(nBoot, bootFun, 1:length(x));
    % isoDO.SalTrendCI(:, i) = prctile(bootTrends(:, 1), [2.5 97.5]);

end

% Display the simplified statistics table with confidence intervals for trends


periods = {'2003-2023'};
time_masks = {isoDO.time > datenum(2003, 1, 1) & isoDO.time < datenum(2024, 1, 1)};

% Initialize arrays to store results
results = struct();

for p = 1:length(periods)
    msk = time_masks{p};
    for i = 1:length(sig_thresh)
        n = sum(~isnan(isoDO.oxygen(i, msk)));
        
        % Calculate trends per decade and uncertainties (remains the same)
        do_trend = isoDO.DOtrend(1, i) * 10;
        do_trend_ci_lower = isoDO.DOtrendCI(1, i) * 10;
        do_trend_ci_upper = isoDO.DOtrendCI(2, i) * 10;
        do_uncertainty =(do_trend_ci_upper - do_trend_ci_lower) / 2;
        
        depth_trend = isoDO.DepthTrend(1, i) * 10;
        depth_trend_ci_lower = isoDO.DepthTrendCI(1, i) * 10;
        depth_trend_ci_upper = isoDO.DepthTrendCI(2, i) * 10;
        depth_uncertainty_trend = (depth_trend_ci_upper - depth_trend_ci_lower) / 2;
        
        temp_trend = isoDO.TempTrend(1, i) * 10;
        temp_trend_ci_lower = isoDO.TempTrendCI(1, i) * 10;
        temp_trend_ci_upper = isoDO.TempTrendCI(2, i) * 10;
        temp_uncertainty_trend = (temp_trend_ci_upper - temp_trend_ci_lower) / 2;
        
        % sal_trend = isoDO.SalTrend(1, i) * 10;
        % sal_trend_ci_lower = isoDO.SalTrendCI(1, i) * 10;
        % sal_trend_ci_upper = isoDO.SalTrendCI(2, i) * 10;
        % sal_uncertainty_trend = (sal_trend_ci_upper - sal_trend_ci_lower) / 2;

        % Calculate mean properties
        mean_depth = mean(isoDO.z(i, msk), 'omitnan');
        mean_oxygen = mean(isoDO.oxygen(i, msk), 'omitnan');
        mean_temp = mean(isoDO.t(i, msk), 'omitnan');
        % mean_sal = mean(isoDO.SA(i, msk), 'omitnan');

        % Calculate uncertainties using bootstrapped confidence intervals
        % Depth
        data_Depth = isoDO.z(i, msk);
        data_Depth = data_Depth(~isnan(data_Depth));
        if length(data_Depth) > 1
            % bootstat_Depth = bootstrp(nBoot, @mean, data_Depth);
            % ci_depth = prctile(bootstat_Depth, [2.5 97.5]);
            depth_uncertainty = nanstd(data_Depth).*2;
        else
            depth_uncertainty = NaN;
        end

        % Dissolved Oxygen
        data_DO = isoDO.oxygen(i, msk);
        data_DO = data_DO(~isnan(data_DO));
        if length(data_DO) > 1
            % bootstat_DO = bootstrp(nBoot, @mean, data_DO);
            % ci_oxygen = prctile(bootstat_DO, [2.5 97.5]);
            % oxygen_uncertainty = (ci_oxygen(2) - ci_oxygen(1)) / 2;
            oxygen_uncertainty = nanstd(data_DO)*2;
        else
            oxygen_uncertainty = NaN;
        end

        % Temperature
        data_Temp = isoDO.t(i, msk);
        data_Temp = data_Temp(~isnan(data_Temp));
        if length(data_Temp) > 1
            % bootstat_Temp = bootstrp(nBoot, @mean, data_Temp);
            % ci_temp = prctile(bootstat_Temp, [2.5 97.5]);
            temp_uncertainty = nanstd(data_Temp)*2;
        else
            temp_uncertainty = NaN;
        end

        % % Salinity
        % data_Sal = isoDO.SA(i, msk);
        % data_Sal = data_Sal(~isnan(data_Sal));
        % if length(data_Sal) > 1
        %     bootstat_Sal = bootstrp(nBoot, @mean, data_Sal);
        %     ci_sal = prctile(bootstat_Sal, [2.5 97.5]);
        %     sal_uncertainty = (ci_sal(2) - ci_sal(1)) / 2;
        % else
        %     sal_uncertainty = NaN;
        % end

        % Store results in the struct
        results(i).sig_thresh = sig_thresh(i);
        results(i).n = n;
        results(i).mean_depth = mean_depth;
        results(i).depth_uncertainty = depth_uncertainty;
        results(i).mean_oxygen = mean_oxygen;
        results(i).oxygen_uncertainty = oxygen_uncertainty;
        results(i).mean_temp = mean_temp;
        results(i).temp_uncertainty = temp_uncertainty;
        % results(i).mean_sal = mean_sal;
        % results(i).sal_uncertainty = sal_uncertainty;
        results(i).do_trend = do_trend;
        results(i).do_uncertainty = do_uncertainty;
        results(i).depth_trend = depth_trend;
        results(i).depth_uncertainty_trend = depth_uncertainty_trend;
        results(i).temp_trend = temp_trend;
        results(i).temp_uncertainty_trend = temp_uncertainty_trend;
        % results(i).sal_trend = sal_trend;
        % results(i).sal_uncertainty_trend = sal_uncertainty_trend;
    end
end

% Print the results
for i = 1:length(results)
    % Format uncertainties for printing
    if isnan(results(i).depth_uncertainty)
        depth_unc_str = 'N/A';
    else
        depth_unc_str = sprintf('± %1.0f', results(i).depth_uncertainty);
    end
    if isnan(results(i).oxygen_uncertainty)
        oxygen_unc_str = 'N/A';
    else
        oxygen_unc_str = sprintf('± %4.1f', results(i).oxygen_uncertainty);
    end
    if isnan(results(i).temp_uncertainty)
        temp_unc_str = 'N/A';
    else
        temp_unc_str = sprintf('± %4.2f', results(i).temp_uncertainty);
    end
    % if isnan(results(i).sal_uncertainty)
    %     sal_unc_str = 'N/A';
    % else
    %     sal_unc_str = sprintf('± %4.2f', results(i).sal_uncertainty);
    % end

    fprintf('%s: %3.2f kg/m^3, n=%4i, mean Depth: %3.0f (%s) m, mean DO: %4.1f (%s) µmol kg^{-1}, mean Temp: %4.2f (%s) °C, DO trend: %4.1f (± %4.1f) µmol kg^{-1} per decade, Depth trend: %3.1f (± %3.1f) m per decade, Temp trend: %4.2f (± %4.2f) °C per decade\n', ...
        periods{p}, results(i).sig_thresh, results(i).n, ...
        results(i).mean_depth, depth_unc_str, ...
        results(i).mean_oxygen, oxygen_unc_str, ...
        results(i).mean_temp, temp_unc_str, ...
        results(i).do_trend, results(i).do_uncertainty, ...
        results(i).depth_trend, results(i).depth_uncertainty_trend, ...
        results(i).temp_trend, results(i).temp_uncertainty_trend)%, ...
        % results(i).sal_trend, results(i).sal_uncertainty_trend);
end


for i=1:length(sig_thresh)
    msk1=find(isoDO.predictY(:,i)<61,1,'first');
    msk2=find(isoDO.predictY_high(:,i)<61,1,'first');
    msk3=find(isoDO.predictY_low(:,i)<61,1,'first');
    fprintf('%4.2f isopycnal hypoxic in %.0f, quickest %.0f, slowest %.0f\n',...
        sig_thresh(i),predictYrs(msk1),predictYrs(msk2),predictYrs(msk3));
end

% save('/Users/samuelstevens/Dropbox/Hakai/gliders/data/isoDO_shelfMn.mat','isoDO');

%% Create glider timeseries by year with uncertainty limits
years = unique(year(allProfs.mtime));  % Extract unique years from mtime
isoO = NaN(length(sig_thresh), length(years));  % Adjust size to match `sig_thresh`
mt = NaN(1, length(years));

% Initialize matrix to store STD
isoO_STD = NaN(length(sig_thresh), length(years));

for i = 1:length(years)-1
    % Mask for each unique year and filter as before
    msk = year(allProfs.mtime) == years(i) & allProfs.missionIdx == 1 & ...
          allProfs.canyonX' > 0 & allProfs.canyonX' < 120;
    mt(i) = mean(allProfs.mtime(msk));

    density_profile = allProfs.s_t(:, msk);
    oxygen_profile = allProfs.oxygen_corrected(:, msk);

    for ii = 1:length(sig_thresh)
        % Pre-allocate an array to store interpolated values
        tmp = NaN(1, sum(msk));

        % Loop through each profile in the mask
        for iii = 1:sum(msk)
            % Extract density and oxygen data for the specific profile
            density_vals = density_profile(:, iii);
            oxygen_vals = oxygen_profile(:, iii);

            % Exclude NaNs from the density and oxygen data
            valid_mask = ~isnan(density_vals) & ~isnan(oxygen_vals);
            density_vals = density_vals(valid_mask);
            oxygen_vals = oxygen_vals(valid_mask);

            % Interpolate if there are enough valid data points
            if numel(density_vals) > 1  % Require at least two points for interpolation
                tmp(iii) = interp1(density_vals, oxygen_vals, sig_thresh(ii));
            end
        end

        % Exclude NaNs from tmp
        tmp_nonan = tmp(~isnan(tmp));

        % Store the mean interpolated oxygen value for the current year and threshold
        isoO(ii, i) = mean(tmp_nonan);  % mean of non-NaN values

        % Calculate STD as uncertainty
        if numel(tmp_nonan) > 0  % Need at least one data point
            % isoO_STD(ii, i) = mean(abs(tmp_nonan - isoO(ii, i)));
            % keyboard
            isoO_STD(ii, i) = nanstd(tmp_nonan);
        else
            isoO_STD(ii, i) = NaN;
        end
    end
end

for i=1:3
    isoDO.residuals(i,:)=isoDO.AmnDO(i,:)-isoDO.fittedy(:, i)';
    isoDO.residualsSD(i,:)=nanstd(isoDO.residuals(i,:));
end


%% How will the frequency of hypoxic days change? 
load S2.mat
load mooringDistsMean.mat

% Parameters
nyrs=90;
T = nyrs * 365; % Total simulation time in days (75 years)
dt = 1; % Time step in days
t = 0:dt:T-1; % Time vector in days

startyr=2016;

% Given model parameters
% O_b0 = 65.8; % Initial bottom oxygen (μmol/kg)
% msk = S2.time > datenum(2017,06,01) & S2.time < datenum(2023,07,01);
% O_b0 = nanmean(S2.CTD280m.oxygen(msk))+((2023-startyr)*0.47);
msk = S2.time > datenum(2016,06,01) & S2.time < datenum(2023,07,01);
O_b0 = nanmean(S2.CTD280m.oxygen(msk));

a = -0.47 / 365; % Daily oxygen decline (μmol/kg per day)
b = 11.9/2; % Interannual variation amplitude (μmol/kg)
T1 = 7 * 365; % Interannual variation period (days)

% Use observed seasonal climatology instead of sinusoidal variation
days_in_year = 365;
clim_time = 1:days_in_year; % Day-of-year indices for climatology
seasonal_cycle = SclimWM(:)-nanmean(SclimWM); % Ensure it is a column vector

% Extend the seasonal cycle to match the full simulation duration
seasonal_variation = repmat(seasonal_cycle, ceil(T/days_in_year), 1);

% Compute oxygen levels over time using the empirical seasonal cycle
% with interannual sinusoid having a phase shift of -pi/2
O_b = O_b0 + a * t + b * sin((2*pi/T1) * t - pi/2) + seasonal_variation';

% Determine hypoxic days (O_b < 61 μmol/kg)
hypoxic_days = O_b < 61;

% Compute number of hypoxic days per year using a loop
Syears = 0:nyrs; % Years from 0 to 75
hypoxic_days_per_year = NaN(size(Syears)); % Initialize array

for i = 1:length(Syears)-1
    start_idx = (i - 1) * 365 + 1; % Start index for the year
    end_idx = i * 365; % End index for the year
    hypoxic_days_per_year(i) = sum(hypoxic_days(start_idx:end_idx));
end

Syrs=year(S2.time);
Syrgrid=2017:2023;
S2_hypoxic_days_per_year=[];
for i=1:length(Syrgrid)
    msk=Syrs==Syrgrid(i);
    S2_hypoxic_days_per_year(i)=...
        sum(S2.CTD280m.oxygen(msk)<=61)/24;
end  

% Plot results
figure;

% Plot oxygen concentration over time
subplot(2,1,1);
plot(startyr + t / 365, O_b, 'b', 'LineWidth', 1.5); % Convert time to years
xlabel('Time (years)');
ylabel('Bottom Oxygen Concentration (\mu mol/kg)');
title('Projected Bottom Oxygen at Scott2 Over 75 Years (Using Seasonal Climatology)');
grid on;
xlim([startyr 2100])

% Plot number of hypoxic days per year
subplot(2,1,2);
bar(startyr + Syears, hypoxic_days_per_year, 'r','facealpha',0.5);
hold on
bar(Syrgrid,S2_hypoxic_days_per_year, 'b','facealpha',0.5);
xlabel('Time (years)');
ylabel('Number of Hypoxic Days per Year');
title('Annual Count of Hypoxic Days (O_b < 61 \mu mol/kg)');
grid on;
xlim([startyr 2100])
legend('model','Scott2');

figure('color','w');
subplot(2,1,1);
plot(startyr + t / 365, O_b, 'b', 'LineWidth', 1.5); % Convert time to years
xlabel('Time (years)');
ylabel('$\mathrm{O_b}$ ($\mu$mol kg$^{-1}$)','interpreter','latex');
title('Projected Bottom Oxygen at Scott2');
grid on;
xlim([startyr 2100])
text(0.01,0.975,'a)','units','normalized');

subplot(2,1,2);
plot(startyr + t / 365, O_b, 'b', 'LineWidth', 1.5); % Convert time to years
hold on
plot(Syrs'+(day(datetime(datestr(S2.time)),'dayofyear')./365)+(hour(S2.time')./24/365),...
    S2.CTD280m.oxygen,'r');
xlim([2016 2024]);
grid on;
legend('Model','Scott2','location','northwest');
ylabel('$\mathrm{O_b}$ ($\mu$mol kg$^{-1}$)','interpreter','latex');
xlabel('Time (years)');
title('Projected Bottom Oxygen at Scott2 vs. Observations');
text(0.01,0.975,'b)','units','normalized');

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/S2model.pdf', ...
    'ContentType', 'vector', 'Resolution', 300);

 %% DO Plotting (Yearly-Averaged Data)
lnes = linspecer(5);lnes(3:4,:)=[]; lnes(1,:)=[0.5 0.5 0.5]; lnes=flipud(lnes);
figure('units', 'centimeters', 'outerposition', [0 0 18.3 9*1.33], 'color', 'w');

Md=load('mooringDistsMean.mat');

ax1=axes('position',[0.05 0.45 0.6 0.79*0.66]);
yrs = 2003:2023;
hold on
for i=1:length(yrs)
    decyrs(i)=datenum(yrs(i),01,01);
end

errorbar(decyrs(15:21), Md.mS, Md.mS_std, 'LineWidth', 0.25, ...
    'CapSize', 2, 'Color', [rgb_x('steel blue') 0.025], 'LineStyle', 'none'); 
plot(decyrs(15:21), Md.mS,'-.', 'Color', rgb_x('steel blue'));
scatter(decyrs(15:21), Md.mS, 40, 'filled', '^', 'MarkerFaceColor', ...
    rgb_x('steel blue'), 'MarkerFaceAlpha', 0.5);

for i = 1:length(sig_thresh)
    % scatter(isoDO.time, isoDO.oxygen(i, :), 2, 'filled', 'MarkerFaceColor', ...
    %     lnes(i, :), 'MarkerFaceAlpha', 0.2);
    plot(decyrs, isoDO.fittedy(:, i), '--', 'color', lnes(i, :));
end

for i = 1:length(sig_thresh)
    jbfill(decyrs(~isnan(isoDO.AmnDO(i, :))), ...
        isoDO.AmnDO(i, ~isnan(isoDO.AmnDO(i, :))) - isoDO.AseDO(i, ~isnan(isoDO.AmnDO(i, :))), ...
        flipud(isoDO.AmnDO(i, ~isnan(isoDO.AmnDO(i, :))) + isoDO.AseDO(i, ~isnan(isoDO.AmnDO(i, :)))), ...
        lnes(i, :), lnes(i, :), 1, 0.085);
    lll(i) = plot(decyrs, isoDO.AmnDO(i, :), '-', 'color', lnes(i, :));
    scatter(decyrs, isoDO.AmnDO(i, :), 20, 'o', 'filled', 'markerfacecolor', ...
        lnes(i, :), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.7);
end

for i = 1:length(sig_thresh)
    errorbar(datenum(2019,01,01):365:datenum(2024,01,01),...
        isoO(i, :), isoO_STD(i,:), 'LineWidth', 0.25, ...
        'CapSize', 2, 'Color', [lnes(i, :) 0.025], 'LineStyle', 'none');
    plot(datenum(2019,01,01):365:datenum(2024,01,01), isoO(i, :),'-.','color',lnes(i, :));
    scatter(datenum(2019,01,01):365:datenum(2024,01,01), isoO(i, :),100,'filled','p', ...
         'MarkerFaceColor', lnes(i, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
end
grid on
axis tight; yl = ylim; ylim([55 195]); axdate(10);
ylabel('O$_2$ ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
yticks([61 75:25:250]);
text(0.025,0.975,'a)','units','normalized');
% xticklabels([]);

ax3=axes('position',[0.05 0.1 0.6 0.275]);
jbfill(decyrs,zeros(size(decyrs))+max(isoDO.residualsSD(i,:)), ...
        zeros(size(decyrs))-max(isoDO.residualsSD(i,:)), ...
        [0.8 0.8 0.8], [0.8 0.8 0.8], 1, 0.4);
hold on
for i=1:3
    errorbar(decyrs,isoDO.residuals(i,:), isoDO.AseDO(i,:), 'LineWidth', 0.25, ...
        'CapSize', 2, 'Color', [lnes(i,:) 0.025], 'LineStyle', 'none');

    plot(decyrs,isoDO.residuals(i,:),'color',lnes(i,:));
    scatter(decyrs,isoDO.residuals(i,:),20,'markerfacecolor',lnes(i,:),...
        'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.7);
end
axdate(10); grid on;
ylabel('O'' ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
xlabel('Year');
axis tight;
text(0.025,0.975,'b)','units','normalized');

% Detrend and calculate interquartile range
for i=1:3
    isoIQR(i)=iqr(isoDO.AmnDO(i,:)-isoDO.fittedy(:,i)');
end

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8);

% Detrend and calculate interquartile range
for i=1:3
    isoIQR(i)=iqr(isoDO.AmnDO(i,:)-isoDO.fittedy(:,i)');
end

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 9);

axes(ax1)
ss1=scatter(NaN, NaN, 20, 'o', 'filled', 'markerfacecolor', ...
        'k', 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.7);
ss2=scatter(NaN, NaN,100,'filled','p', ...
         'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
ss3=scatter(NaN, 40, 'filled', '^', 'MarkerFaceColor', ...
    'k', 'MarkerFaceAlpha', 0.5);
legend([ss1;ss2;ss3],'Ship-based','Glider','Scott2','fontsize',8);

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
%%
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/shelf_isopycnalDO_R1.pdf -dpdf -nofontswap
% exportgraphics(gcf, '/Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/shelf_isopycnalDOV6.pdf', ...
%     'ContentType', 'vector', 'Resolution', 300);
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/poster/ts.png -m6 -nofontswap -transparent

%% Projections
figure('units', 'centimeters', 'outerposition', [0 0 7 13], 'color', 'w');
tiledlayout(2,1,'TileSpacing', 'Compact', 'Padding', 'Compact');

ax2=nexttile;
% ax2=axes('position',[0.75 0.45 0.23 0.79*0.66]);
hold on;
msk=predictYrs<=2100;
for i=1:sum(msk)
    decyrsP(i)=datenum(predictYrs(i),06,01);
end

for i = 1:length(sig_thresh)
    jbfill(decyrsP(msk), ...
        isoDO.predictY_low(msk,i), ...
        flipud(isoDO.predictY_high(msk,i)), ...
        lnes(i, :), lnes(i, :), 1, 0.075);
    plot(decyrsP(msk),isoDO.predictY(msk, i), '--', 'color', lnes(i, :));
end
scatter(decyrsP(1):365*2:max(decyrsP(msk)), ones(size(decyrsP(1):365*2:max(decyrsP(msk)))).*61,...
    3, 's', 'k', 'filled');
grid on;
xticks([datenum(2025,01,01) datenum(2050,01,01) datenum(2075,01,01) datenum(2100,01,01)]);
xticklabels({'2025','2050','2075','2100'});
yticks([61 75:25:250]);
% ax2.YAxis.Visible = 'off';
axis tight;
% xlabel('Year');
ylim([55 195]); 
text(0.025,0.95,'a)','units','normalized');
ylabel('Projected O$_2$ ($\mu$mol kg$^{-1}$)', 'Interpreter', 'latex');
box on 

ax4=nexttile;
bar(startyr + Syears, hypoxic_days_per_year,1, 'k','facealpha',0.3,'edgecolor','none');
hold on
bar(Syrgrid,S2_hypoxic_days_per_year, 1,'facecolor',rgb_x('teal'),...
    'facealpha',0.7,'edgecolor','none');
xlabel('Year');
ylabel('Hypoxic Days per Year');
grid on;
xlim([2010 2100])
text(0.025,0.95,'b)','units','normalized');
scatter(2016,270,30,'s','markerfacecolor','k','markerfacealpha',0.3,...
    'markeredgecolor','none');
text(2020,270,'Model');
scatter(2016,225,30,'s','markerfacecolor',rgb_x('teal'),'markerfacealpha',0.7,...
    'markeredgecolor','none');
text(2020,225,'Scott2');
ylim([0 365]);
yticks(0:120:360)

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
%% 
export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/projections_R1.pdf -dpdf -nofontswap

%% PSD estimates of annual oxygen residuals
figure;

for i = 1:length(sig_thresh)
    % Load residuals data (assuming isoDO.residuals is available)
    residuals = inpaint_nans(isoDO.residuals(i,:)); % Select the isopycnal of interest

    % Define parameters
    N = length(residuals);    % Number of data points
    dt = 1;                   % Sampling interval (years)
    fs = 1/dt;                % Sampling frequency (1/year)

    % Detrend residuals to remove linear trend
    residuals_detrended = detrend(residuals);

    % Compute Power Spectral Density using Welch's method
    window = hamming(floor(N/2));       % Window length: half the series length
    noverlap = floor(length(window)/2); % 50% overlap
    nfft = max(256, 2^nextpow2(N));     % FFT length

    [Pxx, f] = pwelch(residuals_detrended, window, noverlap, nfft, fs);

    % Identify dominant frequency (excluding zero frequency)
    [~, peak_idx] = max(Pxx(2:end));
    dominant_frequency = f(peak_idx + 1);

    % Convert to period in years
    T1_estimated = 1 / dominant_frequency;

    % Plot PSD
    subplot(1,3,i);
    plot(1./f(2:end), Pxx(2:end), 'b', 'LineWidth', 1.5);
    xlabel('Period (Years)');
    ylabel('Power Spectral Density (PSD)');
    title('PSD of Oxygen Residuals');
    grid on;
    xlim([0 20]); % Limit to relevant periods (e.g., 0-20 years)

    % Display estimated dominant period
    fprintf('Estimated Interannual Period (T1): %.2f years\n', T1_estimated);
end

%% simple schematic poster plot

lnes = linspecer(5);lnes(3:4,:)=[]; lnes(1,:)=[0.5 0.5 0.5]; lnes=flipud(lnes);
figure('units', 'centimeters', 'outerposition', [0 0 16 10], 'color', 'w');

ax1=axes('position',[0.1 0.45 0.6 0.79*0.66]);
hold on
Md=load('mooringDistsMean.mat');

for i=1:length(yrs)
    decyrs(i)=datenum(yrs(i),01,01);
end

for i = 3
    % scatter(isoDO.time, isoDO.oxygen(i, :), 2, 'filled', 'MarkerFaceColor', ...
    %     lnes(i, :), 'MarkerFaceAlpha', 0.2);
    plot(decyrs, isoDO.fittedy(:, i), '--', 'color', lnes(i, :),'linewidth',2);
    lll(i) = plot(decyrs, inpaint_nans(isoDO.AmnDO(i, :)), '-', 'color', lnes(i, :),'linewidth',2);
    scatter(decyrs, inpaint_nans(isoDO.AmnDO(i, :)), 40, 'o', 'filled', 'markerfacecolor', ...
        lnes(i, :), 'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.7);
end
line([min(decyrs) max(decyrs)],[61 61],'linestyle','--','color','k')
grid on
axis tight; axdate(4);
ylabel('O$_2$', 'Interpreter', 'latex');
yticks([]);%yticklabels('Hypoxia');
ylim([0 95]); 
yticklabels();

ax2=axes('position',[0.72 0.45 0.23 0.79*0.66]);
hold on;
msk=predictYrs<=2100;
for i=1:sum(msk)
    decyrsP(i)=datenum(predictYrs(i),06,01);
end

for i = 3
    jbfill(decyrsP(msk), ...
        isoDO.predictY_low(msk,i), ...
        flipud(isoDO.predictY_high(msk,i)), ...
        lnes(i, :), lnes(i, :), 1, 0.075);
    plot(decyrsP(msk),isoDO.predictY(msk, i), '--', 'color', lnes(i, :),'linewidth',2);
end
line([min(decyrsP(msk)) max(decyrsP(msk))],[61 61],'linestyle','--','color','k')


grid on;
xticks([datenum(2025,01,01) datenum(2050,01,01) datenum(2075,01,01) datenum(2100,01,01)]);
xticklabels({'2025','2050','2075','2100'});
yticks([61]);yticklabels('Hypoxia');
ax2.YAxis.Visible = 'off';
axis tight;
xlabel('Year');
ylim([0 95]); 

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 20);
set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/shelf_isopycnalDOposter.png -m5 -nofontswap


%% Map
col=[255 214 140]/255; % YELLOW!
lat_lim=[50.5 52.75]; lon_lim=[-131.5 -127];
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

figure('units','centimeters','outerposition',[0 0 15 15],'color','w');
ax1=axes;%('position',[0.1 0.3 0.5 0.5]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
[CS,CH]=m_contour(lon(ilon),lat(ilat),Z',[-100 -500],'linecolor',...
    [0.8 0.8 0.8]);
[CS2,CH2]=m_contour(lon(ilon),lat(ilat),Z',[-200 -200],'linecolor',...
    [0.3 0.3 0.3]);
clabel(CS,CH,'manual','color',[0.8 0.8 0.8],'fontsize',8,...
    'fontname','Latin Modern Roman');
clabel(CS2,CH2,'manual','color',[0.3 0.3 0.3],'fontsize',8,...
    'fontname','Latin Modern Roman');
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
m_grid('box','fancy','tickdir','out','linestyle','none','xtick',3,...
    'YaxisLocation','right');

lnes=linspecer(25);count=0;
for i=1:max(allProfs.missionNum)
    msk=allProfs.missionNum==i;
    if unique(allProfs.missionIdx(msk))==1 && ...
            mean(year(allProfs.mtime(msk)),'omitnan')>=2023 &&...
        mean(month(allProfs.mtime(msk)),'omitnan')>=7 &&...
        mean(month(allProfs.mtime(msk)),'omitnan')<=9
        count=count+1;
        m_plot(allProfs.longitude(msk),allProfs.latitude(msk),...
            'linewidth',0.5,'color',[lnes(count,:) 0.4]);
    end
end
sm(1)=m_scatter(mean(CS09.lon),mean(CS09.lat),125,'filled','p','markerfacecolor',rgb_x('violet'),...
    'markeredgecolor',rgb_x('violet'),'markerfacealpha',0.8);
% sm(2)=m_scatter(isoDO.longitude,isoDO.latitude,20,'filled','markerfacealpha',0.4,...
%     'markerfacecolor',rgb_x('baby blue'));

m_annotation('textarrow',[-129.0631  -129.0631],[51.4114+0.2 51.4114],...
    'string',{'Goose';'Trough'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','center','headlength',5);
 
 sm(2)=m_scatter(-129.47,51.13,75,'^','markeredgecolor','k',...
     'markerfacecolor',rgb_x('steel blue'));
% m_text(-129.47,51.13-0.175,{'Scott2';'(280 m)'},'Horizontalalignment','center');
 sm(3)=m_scatter(-128.24,51.71,75,'^','markeredgecolor','k',...
     'markerfacecolor',rgb_x('rose'));
% m_text(-128.24,51.71+0.2,{'Hak1';'(133 m)'},'Horizontalalignment','center');

legend(sm,'CS09','Scott2','Hak1','location','northwest');

axG=axes('position',[0.1 0.4 0.2 0.2]);
m_proj('lambert','long',[-135 -120],'lat',[48 56]);
[CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -500 -200 0 ],'edgecolor','none');
m_gshhs_i('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
m_grid('linest','none','tickdir','out','xtick',3,'ytick',3);
colormap(m_colmap('blues'));  
clim([-7000 000]);
m_text(-123.3305-1,53.6163, {'\it{British}';'\it{Columbia}'},...
    'horizontalalignment','center');

r=m_rectangle(lon_lim(1),lat_lim(1),...
    diff(lon_lim),diff(lat_lim),0,...
    'edgecolor','k','facecolor',rgb_x('forest green'),'facealpha',0.1,...
    'linewidth',1.5);

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/smallMap.pdf -dpdf -nofontswap

%% CS09 %%%%
load CS09_ctd.mat
sa=gsw_SA_from_SP(CS09.sal,CS09.pr,CS09.lon,CS09.lat);
pt=gsw_pt0_from_t(sa,CS09.temp,CS09.pr);
CS09.AOU=aou(CS09.sal,pt,CS09.ox);

cmonth=month(CS09.mtime);cyr=year(CS09.mtime);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;


figure('units','centimeters','outerposition',[0 0 18 18],'color','w');

lims=std(CS09.ox(:,msk),0,2,'omitnan');
mn=mean(CS09.ox(:,msk),2,'omitnan');
ax1=axes('position',[0.1 0.65 0.2 0.3]);
hold on
z=1:1000;

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(CS09.ox(:,msk),1:1000,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2023);
msk_2023=fliplr(msk_2023);

% Subset the 2023 data
tmpO = CS09.ox(:, msk_2023);
tmpP = CS09.pr(:, msk_2023);

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);
 
% Create mask for data below the threshold
msk2 = tmpO < (mn - 2*lims) & tmpO<61;

% Plot points where msk2 is true in dark red
scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
    'MarkerFaceColor', [0.1 0 0]);
% 
% % Create mask for data in 2022
% msk_2022 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2022);
% msk_2022=fliplr(msk_2022);
% 
% % Subset the 2022 data
% tmpO = CS09.ox(:, msk_2022);
% tmpP = CS09.pr(:, msk_2022);
% 
% % Plot 2022 data in red
% l1(1)=plot(tmpO(:,1), tmpP(:,1),':', 'color', [0.7 0 0], 'linewidth', 1);
% l1(2)=plot(tmpO(:,2), tmpP(:,2),':', 'color', [1 0 0], 'linewidth', 1);
% 
% % Create mask for data below the threshold
% msk2 = tmpO < (mn - 2*lims) & tmpO<61;
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

line([61 61],[0 200],'linestyle','--','color','k');
axis ij
grid on
xlim([40 365]);
xlabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
ylabel('Depth (m)');
legend(l1,datestr(CS09.mtime(msk_2023),1))
text(320,30,'a)');

axzoom=axes('position',[0.1+0.125 0.65+0.05 0.05 0.1]);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;
hold on

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(CS09.ox(:,msk),1:1000,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2023);

% Subset the 2023 data
tmpO = CS09.ox(:, msk_2023);
tmpP = CS09.pr(:, msk_2023);

% Plot 2023 data in red
l(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [1 0 0], 'linewidth', 1.5);
l(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);

% Create mask for data below the threshold
msk2 = tmpO < (mn - 2*lims) & tmpO<61;

% Plot points where msk2 is true in dark red
scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
    'MarkerFaceColor', [0.1 0 0]);

line([61 61],[0 200],'linestyle','--','color','k');
axis ij
grid on
ylim([150 200]);
xlim([50 75]);
box on;


%%%%%%%%%%%%%%%%%% AOU
lims=std(CS09.AOU(:,msk),0,2,'omitnan');
mn=mean(CS09.AOU(:,msk),2,'omitnan');
ax1=axes('position',[0.325 0.65 0.2 0.3]);
hold on
z=1:1000;

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(CS09.AOU(:,msk),1:1000,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2023);
msk_2023=fliplr(msk_2023);

% Subset the 2023 data
tmpO = CS09.AOU(:, msk_2023);
tmpP = CS09.pr(:, msk_2023);

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% Create mask for data below the threshold
msk2 = tmpO > (mn + 2*lims);

% Plot points where msk2 is true in dark red
scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
    'MarkerFaceColor', [0.1 0 0]);
axis ij
grid on
xlim([-80 260]);
xlabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
yticklabels([]);
text(240,10,'b)');

axzoom=axes('position',[0.375 0.65+0.05 0.05 0.1]);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;
hold on

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(CS09.AOU(:,msk),1:1000,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2023);

% Subset the 2023 data
tmpO = CS09.AOU(:, msk_2023);
tmpP = CS09.pr(:, msk_2023);

% Plot 2023 data in red
l(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [1 0 0], 'linewidth', 1.5);
l(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
 
% Create mask for data below the threshold
msk2 = tmpO > (mn + 2*lims);

% Plot points where msk2 is true in dark red
scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
    'MarkerFaceColor', [0.1 0 0]);

axis ij
grid on
ylim([150 200]);
xlim([235 260]);
box on;

%%%%%%%%%%%%%%%%%% s_t
lims=std(CS09.s_t(:,msk),0,2,'omitnan');
mn=mean(CS09.s_t(:,msk),2,'omitnan');
ax1=axes('position',[0.55 0.65 0.2 0.3]);
hold on
z=1:1000;

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(CS09.s_t(:,msk),1:1000,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2023);
msk_2023=fliplr(msk_2023);

% Subset the 2023 data
tmpO = CS09.s_t(:, msk_2023);
tmpP = CS09.pr(:, msk_2023);

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% Create mask for data below the threshold
msk2 = tmpO > (mn + 2*lims);

% Plot points where msk2 is true in dark red
scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
    'MarkerFaceColor', [0.1 0 0]);
axis ij
grid on
xlim([21.9 26.8]);
xlabel('$\sigma_\theta$ (kg m$^{-3}$)','Interpreter','latex');
yticklabels([]);
text(26.6,10,'c)');

axzoom=axes('position',[0.6 0.65+0.05 0.05 0.1]);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;
hold on

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(CS09.s_t(:,msk),1:1000,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 7 & cmonth <= 9 & cyr == 2023);

% Subset the 2023 data
tmpO = CS09.s_t(:, msk_2023);
tmpP = CS09.pr(:, msk_2023);

% Plot 2023 data in red
l(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [1 0 0], 'linewidth', 1.5);
l(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
 
% % Create mask for data below the threshold
% msk2 = tmpO > (mn + 2*lims);
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

axis ij
grid on
ylim([150 200]);
xticks([26.3 26.75]);
set(gca,'XTickLabelRotation',0);
box on;

set(findall(gcf,'-property','fontsize'),'fontsize',12);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/CS09.pdf -dpdf -nofontswap
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/poster/CS09.png -m6 -nofontswap -transparent


%% KC10 %%%%
load KC10.mat
KC=HakaiWaterPropertiesInstrumentP;
KC.time=datenum(1970,01,01,00,00,00)+(KC.time/86400);
KC.pressure=double(KC.pressure);
KC.dissolved_oxygen_ml_l=KC.dissolved_oxygen_ml_l.*43.7;
sa=gsw_SA_from_SP(KC.salinity,KC.pressure,KC.longitude,KC.latitude);
pt=gsw_pt0_from_t(sa,KC.temperature,KC.pressure);
ct=gsw_CT_from_t(sa,KC.temperature,KC.pressure);
KC.s_ti=gsw_sigma0(sa,ct);
KC.AOUi=aou(KC.salinity,pt,KC.dissolved_oxygen_ml_l);

t=unique(KC.time);
KC.ox=NaN(371,length(t));
KC.AOU=NaN(371,length(t));

KC.depth=KC.depth+rand(size(KC.depth))/1e5;
for i=1:length(t)
    msk=KC.time==t(i);
    KC.ox(:,i)=interp1(KC.depth(msk),KC.dissolved_oxygen_ml_l(msk),0:370)';
    KC.AOU(:,i)=interp1(KC.depth(msk),KC.AOUi(msk),0:370)';
    KC.s_t(:,i)=interp1(KC.depth(msk),KC.s_ti(msk),0:370)';

end
t(KC.ox(1,:)>550)=[];
KC.ox(:,KC.ox(10,:)>550)=[];
KC.AOU(:,KC.AOU(10,:)<-550)=[];

cmonth=month(t);cyr=year(t);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;

figure('units','centimeters','outerposition',[0 0 18 18],'color','w');

lims=std(KC.ox(:,msk),0,2,'omitnan');
mn=mean(KC.ox(:,msk),2,'omitnan');
ax1=axes('position',[0.1 0.65 0.2 0.3]);
hold on
z=0:370;

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(KC.ox(:,msk),0:370,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 8 & cmonth <= 10 & cyr == 2015);
msk_2023=fliplr(msk_2023);

% Subset the 2023 data
tmpO = KC.ox(:, msk_2023);
tmpP = repmat([0:370]',1,length(msk_2023));

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.5 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(3)=plot(tmpO(:,3), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% % Create mask for data below the threshold
% msk2 = tmpO < (mn - 2*lims) & tmpO<61;
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

line([61 61],[0 370],'linestyle','--','color','k');
axis ij
grid on
xlim([40 365]);
xlabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
ylabel('Depth (m)');
legend(l1,datestr(t(msk_2023),1))
ylim([0 370]);

axzoom=axes('position',[0.1+0.125 0.65+0.05 0.05 0.1]);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;
hold on

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(KC.ox(:,msk),0:370,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 8 & cmonth <= 10 & cyr == 2015);

% Subset the 2023 data
tmpO = KC.ox  (:, msk_2023);
tmpP = repmat([0:370]',1,length(msk_2023));

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.5 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(3)=plot(tmpO(:,3), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% Create mask for data below the threshold
msk2 = tmpO < (mn - 2*lims) & tmpO<61;

% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

line([61 61],[0 370],'linestyle','--','color','k');
axis ij
grid on
ylim([100 370]);
xlim([55 75]);
box on;


%%%%%%%%%%%%%%%%%% AOU
lims=std(KC.AOU(:,msk),0,2,'omitnan');
mn=mean(KC.AOU(:,msk),2,'omitnan');
ax1=axes('position',[0.325 0.65 0.2 0.3]);
hold on
z=0:370;

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(KC.AOU(:,msk),0:370,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 8 & cmonth <= 10 & cyr == 2015);
msk_2023=fliplr(msk_2023);

% Subset the 2023 data
tmpO = KC.AOU(:, msk_2023);
tmpP = repmat([0:370]',1,length(msk_2023));

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.5 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(3)=plot(tmpO(:,3), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% % Create mask for data below the threshold
% msk2 = tmpO > (mn + 2*lims);
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

axis ij
grid on
xlim([-80 260]);
xlabel('AOU ($\mu$mol kg$^{-1}$)','Interpreter','latex');
yticklabels([]);
ylim([0 370]);

axzoom=axes('position',[0.375 0.65+0.05 0.05 0.1]);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;
hold on

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(KC.AOU(:,msk),0:370,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 8 & cmonth <= 10 & cyr == 2015);

% Subset the 2023 data
tmpO = KC.AOU(:, msk_2023);
tmpP = repmat([0:370]',1,length(msk_2023));

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.5 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(3)=plot(tmpO(:,3), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% % Create mask for data below the threshold
% msk2 = tmpO > (mn + 2*lims);
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

axis ij
grid on
ylim([100 370]);
xlim([190 240]);
box on;

%%%%%%%%%%%%%%%%%% s_t
lims=std(KC.s_t(:,msk),0,2,'omitnan');
mn=mean(KC.s_t(:,msk),2,'omitnan');
ax1=axes('position',[0.55 0.65 0.2 0.3]);
hold on
z=0:370;

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(KC.s_t(:,msk),0:370,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 8 & cmonth <= 10 & cyr == 2015);
msk_2023=fliplr(msk_2023);

% Subset the 2023 data
tmpO = KC.s_t(:, msk_2023);
tmpP = repmat([0:370]',1,length(msk_2023));

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.5 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(3)=plot(tmpO(:,3), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% % Create mask for data below the threshold
% msk2 = tmpO > (mn + 2*lims);
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);
axis ij
grid on
xlim([21.9 26.8]);
xlabel('$\sigma_\theta$ (kg m$^{-3}$)','Interpreter','latex');
yticklabels([]);
ylim([0 370]);

axzoom=axes('position',[0.6 0.65+0.05 0.05 0.1]);
msk=cmonth>=7 & cmonth<=9 & cyr<2023;
hold on

fill([mn(~isnan(mn))-2*lims(~isnan(mn));flipud(mn(~isnan(mn))+2*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.8 0.8 0.8],...
    'edgecolor','none');

fill([mn(~isnan(mn))-1*lims(~isnan(mn));flipud(mn(~isnan(mn))+1*lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.7 0.7 0.7],...
    'edgecolor','none');

plot(KC.s_t(:,msk),0:370,'color',[0.9 0.9 0.9]);

% Create mask for data in 2023
msk_2023 = find(cmonth >= 8 & cmonth <= 10 & cyr == 2015);

% Subset the 2023 data
tmpO = KC.s_t(:, msk_2023);
tmpP = repmat([0:370]',1,length(msk_2023));

% Plot 2023 data in red
l1(1)=plot(tmpO(:,1), tmpP(:,1), 'color', [0.5 0 0], 'linewidth', 1.5);
l1(2)=plot(tmpO(:,2), tmpP(:,2), 'color', [0.7 0 0], 'linewidth', 1.5);
l1(3)=plot(tmpO(:,3), tmpP(:,2), 'color', [1 0 0], 'linewidth', 1.5);

% % Create mask for data below the threshold
% msk2 = tmpO > (mn + 2*lims);
% 
% % Plot points where msk2 is true in dark red
% scatter(tmpO(msk2), tmpP(msk2), 7,'filled', 'color', [0.2 0 0],...
%     'MarkerFaceColor', [0.1 0 0]);

axis ij
grid on
ylim([100 370]);
xlim([25.9 26.35]);

xticks([25.9 26.2]);
set(gca,'XTickLabelRotation',0);
box on;

 
% ax5=axes('position',[0.1 0.3 0.3 0.3]);
% col=[255 214 140]/255; % YELLOW
% lat_lim=[50.5 52]; lon_lim=[-131 -127.5];
% lat_lim2=[51.2 52+0.1];  lon_lim2=[-128.34 -127.65];
% fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
% lat=ncread(fname,'lat');    
% lon=ncread(fname,'lon');
% ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
% ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
% Z=ncread(fname,...
%     'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
%     [ sum(ilon) sum(ilat)],[1 1]);
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
% hold on
% [CS,CH]=m_contour(lon(ilon),lat(ilat),Z',[-200 -200],'linecolor',...
%     [0.3 0.3 0.3]);
% clabel(CS,CH,'labelspacing',600,'color',[0.8 0.8 0.8],'fontsize',6,...
%     'fontname','Latin Modern Roman');
% m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
% m_grid('tickdir','out','linestyle','none','xtick',3,...
%     'YaxisLocation','left');
% m_scatter(mean(KC.longitude,'omitnan'),mean(KC.latitude,'omitnan'),...
%     20,'filled','markerfacecolor','r');

lat_lim2=[51.2 52+0.1];  lon_lim2=[-128.34 -127.65];
ax6=axes('position',[0.275 0.275 0.3 0.3]);
col=[255 214 140]/255; % YELLOW
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim2(1) & lon<=lon_lim2(2);
ilat=lat>=lat_lim2(1) & lat<=lat_lim2(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
m_proj('equidistant','lon',lon_lim2,'lat',lat_lim2);   % Projection
hold on
[CS,CH]=m_contour(lon(ilon),lat(ilat),Z',[-600:100:0],'linecolor',...
    [0.3 0.3 0.3]);
clabel(CS,CH,'labelspacing',600,'color',[0.8 0.8 0.8],'fontsize',6,...
    'fontname','Latin Modern Roman');
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.1);
m_grid('tickdir','out','linestyle','none','xtick',3,...
    'YaxisLocation','left');

m_scatter(mean(KC.longitude,'omitnan'),mean(KC.latitude,'omitnan'),...
    40,'filled','markerfacecolor','r');

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
% export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/figures/KC.pdf -dpdf -nofontswap

%%
% yyaxis right
% 
% for i=1:length(sig_thresh)
%     % jbfill(yrs(~isnan(isoDO.AmnDepth(i,:))),...
%     %     isoDO.AmnDepth(i,~isnan(isoDO.AmnDepth(i,:)))-isoDO.AseDO(i,~isnan(isoDO.AmnDepth(i,:))),...
%     %     flipud(isoDO.AmnDepth(i,~isnan(isoDO.AmnDepth(i,:)))+isoDO.AseDO(i,~isnan(isoDO.AmnDepth(i,:)))),...
%     %     lnes(i,:),lnes(i,:),1,0.1);
%     lll(i)=plot(yrs,isoDO.AmnDepth(i,:),'--','color',lnes(i,:));
%     scatter(yrs,isoDO.AmnDepth(i,:),10,'x','markeredgecolor',...
%         lnes(i,:));
% 
%     % errorbar(yrs,isoDO.AmnDO(i,:),isoDO.AseDO(i,:),'-o',...
%     %     'MarkerSize',5,'markerfacecolor',...
%     %     lnes(i,:),'markeredgecolor','none','color',lnes(i,:),...
%     %     'capsize',4);
% end
% set(gca,'ydir','reverse');
% ylim([0 250])

dd=day(datetime(datestr(isoDO.time)),'dayofyear');
mm=month(isoDO.time);

isoDO.MmnDO=NaN(3,12);
for i=1:12
    msk=mm==i;
    isoDO.MmnDO(:,i)=mean(isoDO.oxygen(:,msk),2,'omitnan');
end
ax2=axes('position',[0.15 0.2 0.6 0.3]); hold on;


for i=1:length(sig_thresh)
    scatter(dd,isoDO.oxygen(i,:),20,'filled','markerfacecolor',...
        lnes(i,:),'markeredgecolor','none','markerfacealpha',0.3);
    plot([1:12].*30,isoDO.MmnDO(i,:),'color',lnes(i,:),'linewidth',1.5);
end

grid on
axis tight; yl=ylim; ylim([55 yl(2)]);
% ylabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
xlabel('Day of Year');

set(findall(gcf,'-property','fontsize'),'fontsize',12);

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


%% OProfiles
% isoDO.Ptime=isoDO.time(isoDO.Poxygen(3,:)>1);
% isoDO.Poxygen(:,isoDO.Poxygen(3,:)<1)=[];
mnth=month(isoDO.time);yr=year(isoDO.time);
msk=mnth>=7 & mnth<=9 & yr<=2022;
% lims=1.96.*std(isoDO.Poxygen(:,msk),0,2,'omitnan')./sqrt(sum(~isnan(isoDO.Poxygen(:,msk))));
lims=1*std(isoDO.Poxygen,0,2,'omitnan');
mn=mean(isoDO.Poxygen(:,msk),2,'omitnan');
z=1:350;

figure('units','centimeters','outerposition',[0 0 10 15],'color','w');
hold on
plot(isoDO.Poxygen(:,msk),1:350,'color',[0.5 0.5 0.5]);
fill([mn(~isnan(mn))-lims(~isnan(mn));flipud(mn(~isnan(mn))+lims(~isnan(mn)))],...
    [z(~isnan(mn)) fliplr(z(~isnan(mn)))],[0.9 0.9 0.9],...
    'edgecolor','none','facealpha',0.4);
% msk=mnth>=7 & mnth<=9 & yr==2022;
% plot(isoDO.Poxygen(:,msk),1:350,'m');
msk=mnth>=7 & mnth<=9 & yr==2023;
plot(isoDO.Poxygen(:,msk),1:350,'r');
axis ij


%% build stats table
clc

msk = isoDO.time > datenum(2002,1,1) & isoDO.time < datenum(2016,1,1) & inShape(AS, isoDO.longitude, isoDO.latitude);
for i = 1:length(sig_thresh)
    n = sum(~isnan(isoDO.oxygen(i, msk)));
    fprintf('%i-2016: %3.2f kg/m3, n=%4i, %3.0f m+-%3.0f,  %4.1f umol/kg+-%4.1f, %2.1f °C+-%2.1f, %4.2f g/kg+-%4.2f\n', ...
        year(min(isoDO.time(msk))), sig_thresh(i), n, ...
        mean(isoDO.z(i, msk), 'omitnan'), 1.96 * std(isoDO.z(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.oxygen(i, msk), 'omitnan'), 1.96 * std(isoDO.oxygen(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.t(i, msk), 'omitnan'), 1.96 * std(isoDO.t(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.SA(i, msk), 'omitnan'), 1.96 * std(isoDO.SA(i, msk), 'omitnan') / sqrt(n));
end

fprintf('\n\n');

msk = isoDO.time > datenum(2016,1,1) & isoDO.time < datenum(2020,1,1) & inShape(AS, isoDO.longitude, isoDO.latitude);
for i = 1:length(sig_thresh)
    n = sum(~isnan(isoDO.oxygen(i, msk)));
    fprintf('2016-2020: %3.2f kg/m3, n=%4i, %3.0f m+-%3.0f,  %4.1f umol/kg+-%4.1f, %2.1f °C+-%2.1f, %4.2f g/kg+-%4.2f\n', ...
        sig_thresh(i), n, ...
        mean(isoDO.z(i, msk), 'omitnan'), 1.96 * std(isoDO.z(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.oxygen(i, msk), 'omitnan'), 1.96 * std(isoDO.oxygen(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.t(i, msk), 'omitnan'), 1.96 * std(isoDO.t(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.SA(i, msk), 'omitnan'), 1.96 * std(isoDO.SA(i, msk), 'omitnan') / sqrt(n));
end

fprintf('\n\n');

msk = isoDO.time > datenum(2020,1,1) & inShape(AS, isoDO.longitude, isoDO.latitude);
for i = 1:length(sig_thresh)
    n = sum(~isnan(isoDO.oxygen(i, msk)));
    fprintf('2020 onward: %3.2f kg/m3, n=%4i, %3.0f m+-%3.0f,  %4.1f umol/kg+-%4.1f, %2.1f °C+-%2.1f, %4.2f g/kg+-%4.2f\n', ...
        sig_thresh(i), n, ...
        mean(isoDO.z(i, msk), 'omitnan'), 1.96 * std(isoDO.z(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.oxygen(i, msk), 'omitnan'), 1.96 * std(isoDO.oxygen(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.t(i, msk), 'omitnan'), 1.96 * std(isoDO.t(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.SA(i, msk), 'omitnan'), 1.96 * std(isoDO.SA(i, msk), 'omitnan') / sqrt(n));
end


%% lowest n-year period
yrs=year(isoDO.time);
wn=6;

triad=[];
for i=2002:2023-wn
    msk=yrs>=i & yrs<i+wn;
    triad(:,i-2001)=mean(isoDO.oxygen(:,msk),2,'omitnan');
end
x=2002+wn/2:1:2023-wn/2;

figure;
plot(x,triad(1,:)); hold on
plot(x,triad(2,:));
plot(x,triad(3,:));


%% Detrend and calculate interquartile range
for i=1:3
    isoIQR(i)=iqr(isoDO.AmnDO(i,:)-isoDO.fittedy(:,i)');
end


%%
clc
msk = isoDO.time > datenum(2002,1,1) & isoDO.time < datenum(2017,1,1) & inShape(AS, isoDO.longitude, isoDO.latitude);
for i = 1:length(sig_thresh)
    n = sum(~isnan(isoDO.oxygen(i, msk)));
    fprintf('until 2017: %3.2f kg/m3, n=%4i, %3.0f m+-%3.0f,  %4.1f umol/kg+-%4.1f, %2.1f °C+-%2.1f, %4.2f g/kg+-%4.2f\n', ...
        sig_thresh(i), n, ...
        mean(isoDO.z(i, msk), 'omitnan'), 1.96 * std(isoDO.z(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.oxygen(i, msk), 'omitnan'), 1.96 * std(isoDO.oxygen(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.t(i, msk), 'omitnan'), 1.96 * std(isoDO.t(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.SA(i, msk), 'omitnan'), 1.96 * std(isoDO.SA(i, msk), 'omitnan') / sqrt(n));
end

fprintf('\n\n');

msk = isoDO.time > datenum(2017,1,1) & inShape(AS, isoDO.longitude, isoDO.latitude);
for i = 1:length(sig_thresh)
    n = sum(~isnan(isoDO.oxygen(i, msk)));
    fprintf('2017-2023: %3.2f kg/m3, n=%4i, %3.0f m+-%3.0f,  %4.1f umol/kg+-%4.1f, %2.1f °C+-%2.1f, %4.2f g/kg+-%4.2f\n', ...
        sig_thresh(i), n, ...
        mean(isoDO.z(i, msk), 'omitnan'), 1.96 * std(isoDO.z(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.oxygen(i, msk), 'omitnan'), 1.96 * std(isoDO.oxygen(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.t(i, msk), 'omitnan'), 1.96 * std(isoDO.t(i, msk), 'omitnan') / sqrt(n), ...
        mean(isoDO.SA(i, msk), 'omitnan'), 1.96 * std(isoDO.SA(i, msk), 'omitnan') / sqrt(n));
end

%%
 
% %% How will the frequency of hypoxic days chamnge? 
% 
% load S2.mat
% load mooringDistsMean.mat
% 
% % Parameters
% T = 75 * 365; % Total simulation time in days (75 years)
% dt = 1; % Time step in days
% t = 0:dt:T; % Time vector in days
% 
% % Given model parameters
% % O_b0 = 65.8; % Initial bottom oxygen (μmol/kg)
% msk=S2.time>datenum(2017,06,01) & S2.time<datenum(2023,07,01);
% O_b0 = nanmean(S2.CTD280m.oxygen(msk));
% a = -0.47 / 365; % Daily oxygen decline (μmol/kg per day)
% b = 11.9; % Interannual variation amplitude (μmol/kg)
% T1 = 7 * 365; % Interannual variation period (days)
% c = 17.0; % Seasonal variation amplitude (μmol/kg)
% T2 = 365; % Seasonal variation period (days)
% 
% 
% % Compute oxygen levels over time
% O_b = O_b0 + a * t + b * sin((2*pi/T1) * t) + c * sin((2*pi/T2) * t);
% 
% % Determine hypoxic days (O_b < 61 μmol/kg)
% hypoxic_days = O_b < 61;
% 
% % Compute number of hypoxic days per year using a loop
% years = 0:75; % Years from 0 to 75
% hypoxic_days_per_year = NaN(size(years)); % Initialize array
% 
% for i = 1:length(years)-1
%     start_idx = (i - 1) * 365 + 1; % Start index for the year
%     end_idx = i * 365; % End index for the year
%     hypoxic_days_per_year(i) = sum(hypoxic_days(start_idx:end_idx));
% end
% 
% % Plot results
% figure;
% 
% % Plot oxygen concentration over time
% subplot(2,1,1);
% plot(2025+t / 365, O_b, 'b', 'LineWidth', 1.5); % Convert time to years
% xlabel('Time (years)');
% ylabel('Bottom Oxygen Concentration (\mu mol/kg)');
% title('Projected Bottom Oxygen at Scott2 Over 75 Years (Daily Resolution)');
% grid on;
% xlim([2025 2100])
% 
% % Plot number of hypoxic days per year
% subplot(2,1,2);
% bar(2025+years, hypoxic_days_per_year, 'r');
% xlabel('Time (years)');
% ylabel('Number of Hypoxic Days per Year');
% title('Annual Count of Hypoxic Days (O_b < 61 \mu mol/kg)');
% grid on;
% xlim([2025 2100])