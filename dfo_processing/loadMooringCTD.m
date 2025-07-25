clear

% Specify the directory containing the NetCDF files
dataDirectory = '/Users/samuelstevens/Dropbox/Hakai/moorings/data/scott2/';
clc
% Get a list of all files in the directory
files = dir(fullfile(dataDirectory, '*ctd*.nc'));

% Check if there are any files to process
if isempty(files)
    error('No files found with the phrase "ctd" in the filename in the specified directory.');
end

% Initialize a structure array to hold all data
allData = struct();

% Loop through each file and read the data
for k = 1:length(files)
    filename = fullfile(dataDirectory, files(k).name);
    fprintf('Reading file: %s\n', filename);
    
    % Use the readMooringCTDData function to read data from the current file
    try
        data = readMooringCTDData(filename);
        allData.(files(k).name(1:33)) = data;
    catch ME
        warning('Failed to read file %s due to the error: %s', filename, ME.message);
        continue;
    end
end

% Now allData contains the data from each file, accessible by the file name
disp('All files have been processed.');

%% Interpolation of all data to a common time grid

flds=fieldnames(allData);
mt200m=[]; mt280m=[];
s200m=[]; s280m=[];
t200m=[]; t280m=[];
o200m=[]; o280m=[];
s_t200m=[]; s_t280m=[];
aou280m=[];
spice280m=[];


l=[]; ll=[];

for i = 1:length(flds)

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-200)<10
        mt200m=[mt200m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t200m=[t200m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s200m=[s200m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t200m=[s_t200m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o200m=[o200m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 

        l=[l;allData.(flds{i}).longitude];
        ll=[ll;allData.(flds{i}).latitude];

    end

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-280)<10
        mt280m=[mt280m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t280m=[t280m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s280m=[s280m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t280m=[s_t280m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o280m=[o280m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 

        aou280m=[aou280m;reshape(allData.(flds{i}).AOU',...
            [length(allData.(flds{i}).time') 1])]; 

        spice280m=[spice280m;reshape(allData.(flds{i}).spice',...
            [length(allData.(flds{i}).time') 1])]; 
    end
end

mtGrid=datenum(2017,06,01):1/24:today';  % hourly intervals
S2.CTD200m.temperature=interp1(mt200m,t200m,mtGrid);
S2.CTD280m.temperature=interp1(mt280m,t280m,mtGrid);

S2.CTD200m.salinity = interp1(mt200m, s200m, mtGrid);
S2.CTD280m.salinity = interp1(mt280m, s280m, mtGrid);

S2.CTD200m.sigma_theta = interp1(mt200m, s_t200m, mtGrid);
S2.CTD280m.sigma_theta = interp1(mt280m, s_t280m, mtGrid);

S2.CTD200m.oxygen = interp1(mt200m, o200m, mtGrid);
S2.CTD280m.oxygen = interp1(mt280m, o280m, mtGrid);

S2.CTD280m.aou = interp1(mt280m, aou280m, mtGrid);

S2.CTD280m.spice = interp1(mt280m, spice280m, mtGrid);

S2.time=mtGrid;
S2.longitude=double(mean(l,'omitnan'));
S2.latitude=double(mean(ll,'omitnan'));

save('/Users/samuelstevens/Dropbox/Hakai/moorings/data/S2.mat',"S2");
clearvars S2

mtGrid=datenum(2016,06,01):1/24:datenum(2017,06,01)-1;  % hourly intervals
S2.CTD200m.temperature=interp1(mt200m,t200m,mtGrid);
S2.CTD280m.temperature=interp1(mt280m,t280m,mtGrid);

S2.CTD200m.salinity = interp1(mt200m, s200m, mtGrid);
S2.CTD280m.salinity = interp1(mt280m, s280m, mtGrid);

S2.CTD200m.sigma_theta = interp1(mt200m, s_t200m, mtGrid);
S2.CTD280m.sigma_theta = interp1(mt280m, s_t280m, mtGrid);

S2.CTD200m.oxygen = interp1(mt200m, o200m, mtGrid);
S2.CTD280m.oxygen = interp1(mt280m, o280m, mtGrid);

S2.CTD280m.aou = interp1(mt280m, aou280m, mtGrid);

S2.CTD280m.spice = interp1(mt280m, spice280m, mtGrid);

S2.time=mtGrid;
S2.longitude=double(mean(l,'omitnan'));
S2.latitude=double(mean(ll,'omitnan'));

save('/Users/samuelstevens/Dropbox/Hakai/moorings/data/S2_17.mat',"S2");
clear
%%
% Specify the directory containing the NetCDF files
dataDirectory = '/Users/samuelstevens/Dropbox/Hakai/moorings/data/scott3/';
clc
% Get a list of all files in the directory
files = dir(fullfile(dataDirectory, '*ctd*.nc'));

% Check if there are any files to process
if isempty(files)
    error('No files found with the phrase "ctd" in the filename in the specified directory.');
end

% Initialize a structure array to hold all data
allData = struct();

% Loop through each file and read the data
for k = 1:length(files)
    filename = fullfile(dataDirectory, files(k).name);
    fprintf('Reading file: %s\n', filename);
    
    % Use the readMooringCTDData function to read data from the current file
    try
        data = readMooringCTDData(filename);
        allData.(files(k).name(1:33)) = data;
    catch ME
        warning('Failed to read file %s due to the error: %s', filename, ME.message);
        continue;
    end
end

% Now allData contains the data from each file, accessible by the file name
disp('All files have been processed.');

%% Interpolation of all data to a common time grid

flds=fieldnames(allData);
mt150m=[]; mt225m=[];
s150m=[]; s225m=[];
t150m=[]; t225m=[];
o150m=[]; o225m=[];
s_t150m=[]; s_t225m=[];
l=[]; ll=[];

for i = 1:length(flds)

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-150)<10
        mt150m=[mt150m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t150m=[t150m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s150m=[s150m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t150m=[s_t150m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])];
        

        o150m=[o150m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 
        l=[l;allData.(flds{i}).longitude];
        ll=[ll;allData.(flds{i}).latitude];

    end

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-225)<10
        mt225m=[mt225m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t225m=[t225m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s225m=[s225m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t225m=[s_t225m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o225m=[o225m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 
    end
end

mtGrid=datenum(2017,06,01):1/24:today';  % hourly intervals
S3.CTD150m.temperature=interp1(mt150m,t150m,mtGrid);
S3.CTD225m.temperature=interp1(mt225m,t225m,mtGrid);

S3.CTD150m.salinity = interp1(mt150m, s150m, mtGrid);
S3.CTD225m.salinity = interp1(mt225m, s225m, mtGrid);

S3.CTD150m.sigma_theta = interp1(mt150m, s_t150m, mtGrid);
S3.CTD225m.sigma_theta = interp1(mt225m, s_t225m, mtGrid);

S3.CTD150m.oxygen = interp1(mt150m, o150m, mtGrid);
S3.CTD225m.oxygen = interp1(mt225m, o225m, mtGrid);

S3.time=mtGrid;
S3.longitude=double(mean(l,'omitnan'));
S3.latitude=double(mean(ll,'omitnan'));


save('/Users/samuelstevens/Dropbox/Hakai/moorings/data/S3.mat',"S3");
clear

%%
% Specify the directory containing the NetCDF files
dataDirectory = '/Users/samuelstevens/Dropbox/Hakai/moorings/data/scott1/';
clc
% Get a list of all files in the directory
files = dir(fullfile(dataDirectory, '*ctd*.nc'));

% Check if there are any files to process
if isempty(files)
    error('No files found with the phrase "ctd" in the filename in the specified directory.');
end

% Initialize a structure array to hold all data
allData = struct();

% Loop through each file and read the data
for k = 1:length(files)
    filename = fullfile(dataDirectory, files(k).name);
    fprintf('Reading file: %s\n', filename);
    
    % Use the readMooringCTDData function to read data from the current file
    try
        data = readMooringCTDData(filename);
        allData.(files(k).name(1:33)) = data;
    catch ME
        warning('Failed to read file %s due to the error: %s', filename, ME.message);
        continue;
    end
end

% Now allData contains the data from each file, accessible by the file name
disp('All files have been processed.');

%% Interpolation of all data to a common time grid

flds=fieldnames(allData);
mt103m=[]; mt015m=[];
s103m=[]; s015m=[];
t103m=[]; t015m=[];
o103m=[]; o015m=[];
s_t103m=[]; s_t015m=[];
l=[]; ll=[];

for i = 1:length(flds)

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-103)<10
        mt103m=[mt103m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t103m=[t103m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s103m=[s103m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t103m=[s_t103m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o103m=[o103m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 

        l=[l;allData.(flds{i}).longitude];
        ll=[ll;allData.(flds{i}).latitude];
    end

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-015)<10
        mt015m=[mt015m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t015m=[t015m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s015m=[s015m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t015m=[s_t015m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o015m=[o015m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 
    end
end

mtGrid=datenum(2017,06,01):1/24:today';  % hourly intervals
S1.CTD103m.temperature=interp1(mt103m,t103m,mtGrid);
S1.CTD015m.temperature=interp1(mt015m,t015m,mtGrid);

S1.CTD103m.salinity = interp1(mt103m, s103m, mtGrid);
S1.CTD015m.salinity = interp1(mt015m, s015m, mtGrid);

S1.CTD103m.sigma_theta = interp1(mt103m, s_t103m, mtGrid);
S1.CTD015m.sigma_theta = interp1(mt015m, s_t015m, mtGrid);

S1.CTD103m.oxygen = interp1(mt103m, o103m, mtGrid);
S1.CTD015m.oxygen = interp1(mt015m, o015m, mtGrid);

S1.time=mtGrid;
S1.longitude=double(mean(l,'omitnan'));
S1.latitude=double(mean(ll,'omitnan'));


save('/Users/samuelstevens/Dropbox/Hakai/moorings/data/S1.mat',"S1");


%%
% Specify the directory containing the NetCDF files
dataDirectory = '/Users/samuelstevens/Dropbox/Hakai/moorings/data/hak1/';
clc
% Get a list of all files in the directory
files = dir(fullfile(dataDirectory, '*ctd*.nc'));

% Check if there are any files to process
if isempty(files)
    error('No files found with the phrase "ctd" in the filename in the specified directory.');
end

% Initialize a structure array to hold all data
allData = struct();

% Loop through each file and read the data
for k = 1:length(files)
    filename = fullfile(dataDirectory, files(k).name);
    fprintf('Reading file: %s\n', filename);
    
    % Use the readMooringCTDData function to read data from the current file
    try
        data = readMooringCTDData(filename);
        allData.(files(k).name(1:31)) = data;
    catch ME
        % keyboard
        warning('Failed to read file %s due to the error: %s', filename, ME.message);
        continue;
    end
end

% Now allData contains the data from each file, accessible by the file name
disp('All files have been processed.');
keyboard
%% Interpolation of all data to a common time grid

flds=fieldnames(allData);
mt133m=[]; mt075m=[];
s133m=[]; s075m=[];
t133m=[]; t075m=[];
o133m=[]; o075m=[];
s_t133m=[]; s_t075m=[];
aou133m=[];
spice133m=[];

l=[]; ll=[];

for i = 1:length(flds)

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-133)<10
        mt133m=[mt133m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t133m=[t133m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s133m=[s133m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t133m=[s_t133m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o133m=[o133m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])];
        
        aou133m=[aou133m;reshape(allData.(flds{i}).AOU',...
            [length(allData.(flds{i}).time') 1])];

        spice133m=[spice133m;reshape(allData.(flds{i}).spice',...
            [length(allData.(flds{i}).time') 1])];

        l=[l;allData.(flds{i}).longitude];
        ll=[ll;allData.(flds{i}).latitude];
    end

    if abs(mean(allData.(flds{i}).pressure,'omitnan')-075)<10
        mt075m=[mt075m;reshape(allData.(flds{i}).time',...
            [length(allData.(flds{i}).time') 1])];
        
        t075m=[t075m;reshape(allData.(flds{i}).temperature',...
            [length(allData.(flds{i}).time') 1])]; 

        s075m=[s075m;reshape(allData.(flds{i}).SA',...
            [length(allData.(flds{i}).time') 1])]; 

        s_t075m=[s_t075m;reshape(allData.(flds{i}).sigma_theta',...
            [length(allData.(flds{i}).time') 1])]; 

        o075m=[o075m;reshape(allData.(flds{i}).oxygen',...
            [length(allData.(flds{i}).time') 1])]; 
    end
end

mtGrid=datenum(2017,06,01):1/24:today';  % hourly intervals
H1.CTD133m.temperature=interp1(mt133m,t133m,mtGrid);
H1.CTD075m.temperature=interp1(mt075m,t075m,mtGrid);

H1.CTD133m.salinity = interp1(mt133m, s133m, mtGrid);
H1.CTD075m.salinity = interp1(mt075m, s075m, mtGrid);

H1.CTD133m.sigma_theta = interp1(mt133m, s_t133m, mtGrid);
H1.CTD075m.sigma_theta = interp1(mt075m, s_t075m, mtGrid);

H1.CTD133m.oxygen = interp1(mt133m, o133m, mtGrid);
H1.CTD075m.oxygen = interp1(mt075m, o075m, mtGrid);

H1.CTD133m.aou = interp1(mt133m, aou133m, mtGrid);

H1.CTD133m.spice = interp1(mt133m, spice133m, mtGrid);

H1.time=mtGrid;
H1.longitude=double(mean(l,'omitnan'));
H1.latitude=double(mean(ll,'omitnan'));

H1.CTD133m.oxygen(H1.time>738180 & H1.time<738312)=NaN;

save('/Users/samuelstevens/Dropbox/Hakai/moorings/data/H1.mat',"H1");