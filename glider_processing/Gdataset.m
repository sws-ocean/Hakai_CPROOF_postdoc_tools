%% Building a glider dataset

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));
%%
clear

%% Read in and store data

% Initialize a structure to store the variables
flist = dir('/Users/samst/Dropbox/Hakai/gliders/data/grids/*grid_delayed.nc');

varList = {'depth', 'profile', 'time', 'longitude', 'latitude',...
    'profile_time_start', 'profile_time_end', 'heading', 'pitch',...
    'roll', 'conductivity', 'temperature', 'pressure', 'chlorophyll',...
    'cdom', 'backscatter_700', 'oxygen_concentration',...
    'temperature_oxygen', 'distance_over_ground', 'profile_index',...
    'profile_direction', 'salinity', 'potential_density', 'density',...
    'potential_temperature', 'waypoint_latitude', 'waypoint_longitude',...
    'oxygen_saturation', 'u', 'v', 'temperature_optics'};

% Initialize all possible variables in the structure with empty arrays
allData = struct();
for i = 1:length(varList)
    allData.(varList{i}) = {};
end
clc
count=0;

% Loop through each file in the file list
for k = 1:length(flist)
    fileList{k} = [flist(k).folder, '/', flist(k).name];
    fileInfo = ncinfo(fileList{k});

    % Determine the maximum array size for initializing NaN arrays
    % (Assuming first file has all representative sizes; adjust as needed)
    exampleFile = ncread(fileList{k}, varList{12});
    defaultSize = size(exampleFile); % adjust based on your data structure needs

    % Initialize NaN arrays for each variable in varList if file does not contain them
    presentVars = {fileInfo.Variables.Name};

    for j = 1:length(varList)
        if ismember(varList{j}, presentVars)
            varData = ncread(fileList{k}, varList{j});
        else
            varData = NaN(defaultSize);
        end
        if isempty(allData.(varList{j}))
            allData.(varList{j}) = {varData};
        else
            allData.(varList{j}){end+1} = varData;
        end
    end

    allData.optode{1,k}=ncreadatt(fileList{k}, '/', 'oxygen');

    % Matlab time
    if double(allData.time{1,k}(1))<1
        
        baseDate=fileInfo.Variables(3).Attributes(6).Value(19:40);
        % Convert base date to datetime
        baseDatetime = datetime(baseDate, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SS', 'TimeZone', 'UTC');

        % Calculate dates by adding nanoseconds converted to days
        dates = baseDatetime + seconds(double(allData.time{1,k}) / 1e9);

        % Convert datetime array to MATLAB datenum
        allData.mtime{1,k} = datenum(dates);

        % allData.mtime{1,k}=double(allData.time{1,k})/1e9+datenum(...
        %     fileInfo.Variables(3).Attributes(6).Value(19:40));
        % count=count+1;
        % nanolist{count}=fileList{k};
        % fprintf('%s \n',fileInfo.Variables(3).Attributes(6).Value)
        % figure;
        % histogram(allData.mtime{1,k});
    else
        allData.mtime{1,k}=double(allData.time{1,k})/86400+datenum('1970', 'yyyy');
        matlab_time = datetime(double(allData.time{1,k}),...
            'ConvertFrom', 'posix', 'TimeZone', 'UTC');
        allData.mtime{1,k}=datenum(matlab_time);
    end
end

% Build a canyon axis
lat_lim=[50.5 52]; lon_lim=[-130.5 -127];
[z.ELEV,z.LON,z.LAT]=m_etopo2([lon_lim lat_lim]);

tempAX=[-128.269656332815	-128.297705598052	-128.318742546980	-128.360816444836	-128.395878026382	-128.434445766083	-128.494050454712	-128.546642827031	-128.641309097206	-128.750000000000	-128.841160112020	-129.026986494216	-129.174245136711	-129.360071518906	-129.658094962050    -129.8369  -130.1420 ;...
    51.7020848573519	51.6713606437454	51.6362472567667	51.5923555230432	51.5528529626920	51.5133504023409	51.4848207754206	51.4475128017557	51.4123994147769	51.4102048280907	51.4058156547184	51.3443672275055	51.3092538405267	51.2543891733724	51.1556327724945    51.0898 50.9880];
tempX=[0 cumsum(gsw_distance(tempAX(1,:),tempAX(2,:)))]/1000; % km

canyonX=0:0.1:max(tempX);
canyonAX=[];
canyonAX(:,1)=interp1(tempX,tempAX(1,:),canyonX);
canyonAX(:,2)=interp1(tempX,tempAX(2,:),canyonX);

canyonZ=interp2(z.LON,z.LAT,z.ELEV,canyonAX(:,1),canyonAX(:,2));

%% Add variables to data structure and remove weird velocities

lat_lim=[50.5 52]; lon_lim=[-130.5 -127];
col=[255 214 140]/255; % YELLOW!
clc

% Loop through each file in the file list
for k = 1:length(flist)

    % Find the closest axis point for every lat/lon point
    for i=1:length(allData.longitude{1,k})

        dd=gsw_distance([ones(length(canyonX),1)*allData.longitude{1,k}(i) ...
            canyonAX(:,1)],...
            [ones(length(canyonX),1)*allData.latitude{1,k}(i) ...
            canyonAX(:,2)]);

        [allData.canyonY{1,k}(i),allData.canyonIDX{1,k}(i)]=min(dd);

        % Add transect direction: 1=offshore, 2=inshore, 3=other
        tmpH=median(allData.heading{1,k}(i,:),'omitnan');
        if tmpH < 90
            allData.transectDir{1,k}(i)=2;
        elseif tmpH > 180 && tmpH < 270
            allData.transectDir{1,k}(i)=1;
        else
            allData.transectDir{1,k}(i)=3;
        end
    end

    % remove odd velocites
    if size(allData.u{1,k},2)~=length(allData.depth{1,k})
        allData.u{1,k}=NaN(size(allData.temperature{1,k}));
        allData.v{1,k}=NaN(size(allData.temperature{1,k}));
    end

    allData.missionNum{1,k}=ones(1,size(allData.temperature{1,k},1)).*k;
    allData.mission{1,k} = repmat({string(flist(k).name)},...
        size(allData.temperature{1, k}, 1), 1);
     allData.optode{1,k} = repmat({string(allData.optode{1,k})},...
        size(allData.temperature{1, k}, 1), 1);

    % % add optode info
    % msk=strcmp(vertcat(fileInfo.Attributes.Name),)
    % 
    % Create index of Calvert Line transects: 0=no, 1=yes, 2=kinda
    figure;
    m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
    m_etopo2('contourf',[-5000 -300:20:0],'edgecolor','none');
    hold on
    m_gshhs_i('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
    clim([-500 0]);
    cm=(m_colmap('blues',34));cm([1 end-1:end],:)=[];
    colormap(cm);
    m_plot(allData.longitude{1,k},allData.latitude{1,k},'linewidth',0.5,...
        'color','r','linewidth',2);

    allData.missionIdx{1,k}=ones(1,size(allData.temperature{1,k},1)).*...
        input('Is this mission through the canyon?\n 0=no, 1=yes, 2=maybe: ');
    close all

end


%% Create allProfs.mat file

% Initialize the output structure
allProfs = struct();

varList=fieldnames(allData);
varList(ismember(varList, {'time';'profile';'profile_time_start';...
    'profile_time_end';'profile_index';'waypoint_latitude';...
    'waypoint_longitude'}))=[];

% Loop through the variables and process
for i = 1:length(varList)
    tmp=allData.(varList{i});
    if size(allData.(varList{i}){1,k},1)>1
        % keyboard
        tmp=vertcat(tmp{:});
        allProfs.(varList{i})=tmp';
    else
        tmp=horzcat(tmp{:});
        allProfs.(varList{i})=tmp;
    end
end

% add correctly sized Z
allProfs.z=ones(size(allProfs.temperature)).*allData.depth{1,1};
allProfs.depth=[];

% TEOS-10 variables
allProfs.SA=gsw_SA_from_SP(allProfs.salinity,...
    allProfs.pressure,allProfs.longitude,allProfs.latitude);
allProfs.CT=gsw_CT_from_t(allProfs.SA,allProfs.temperature,...
    allProfs.pressure);
allProfs.s_t=gsw_sigma0(allProfs.SA,allProfs.CT);
[allProfs.N2,allProfs.midZ]=gsw_Nsquared(allProfs.SA,allProfs.CT,...
    allProfs.pressure,allProfs.latitude);

allProfs=rmfield(allProfs,'depth');

% Sort all fields by mtime
[~, sortedIndices] = sort(allProfs.mtime);

varList = fieldnames(allProfs); % Get list of all fields in allProfs
for i = 1:length(varList)
    allProfs.(varList{i}) = allProfs.(varList{i})(:, sortedIndices);
end

% Change mission number to be time sorted
tmpM=unique(allProfs.missionNum,'stable');
for i=1:length(tmpM)
    msk=allProfs.missionNum==tmpM(i);
    allProfs.missionNum(msk)=i+max(tmpM);
end
allProfs.missionNum=allProfs.missionNum-max(tmpM);


%%
canyonX=canyonX';
allProfs.canyonX=canyonX(allProfs.canyonIDX);
save("/Users/samuelstevens/Dropbox/Hakai/gliders/data/glider_allProfs_May25.mat",'allProfs','canyonAX','canyonX',...
    'canyonZ','-v7.3');
