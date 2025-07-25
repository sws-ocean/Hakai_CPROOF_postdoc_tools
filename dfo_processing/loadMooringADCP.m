clear

% Specify the directory containing the NetCDF files
dataDirectory = '/Users/samuelstevens/Dropbox/Hakai/moorings/data/scott2/';
clc
% Get a list of all files in the directory with the pattern *adcp.nc
files = dir(fullfile(dataDirectory, '*adcp.nc'));

% Check if there are any files to process
if isempty(files)
    error('No files found with the phrase "adcp" in the filename in the specified directory.');
end

% Initialize a structure array to hold all data
allData = struct();

% Loop through each file and read the data
for k = 1:length(files)
    filename = fullfile(dataDirectory, files(k).name);
    fprintf('Reading file: %s\n', filename);
    
    % Use the updated readMooringADCPData function to read data from the current file
    try
        data = readMooringADCPData(filename);
        % Using the first 33 characters of the file name as the structure field name
        % Ensure the field name does not contain invalid characters
        fieldName = matlab.lang.makeValidName(files(k).name(1:33));
        allData.(fieldName) = data;
    catch ME
        warning('Failed to read file %s due to the error: %s', filename, ME.message);
        continue;
    end
end

% Now allData contains the data from each file, accessible by the file name
disp('All files have been processed.');


%% Calculate backscatter
% Constants
kc = 0.5; % Scale factor for converting RSSI counts to dB - adjust as needed
a = 24.5/1000; % Attenuation coefficient in dB/m - adjust based on your specific environment or ADCP documentation

% Loop through each set of data in allData
fields = fieldnames(allData);
for i = 1:length(fields)

    % Current dataset
    currentData = allData.(fields{i});

    R0 = currentData.distance(1); % Reference distance from the transducer in meters - adjust as needed


    % Extract echo intensity from four beams
    E1 = currentData.TNIHCE01;
    E2 = currentData.TNIHCE02;
    E3 = currentData.TNIHCE03;
    E4 = currentData.TNIHCE04;
    
    % Calculate the average echo intensity across all beams
    E = nanmean(cat(3, E1, E2, E3, E4), 3);
    % 
    % % Assuming distance and time dimensions match those of E
    % R = repmat(currentData.distance, [size(E, 1), 1]);  % Replicate distance to match the time dimension of E
    % 
    % % Assume a simplistic noise level estimation or retrieve it if it's calculated elsewhere
    % Enoise = mean(E, 'omitnan');  % This is a simplistic estimation, consider more specific noise level data
    % 
    % % Calculate corrected echo intensity in dB
    % EdB = kc * (E - Enoise);
    % 
    % % Calculate backscatter strength using equation (3)
    % sy = EdB - repmat(EdB(:, 1), [1, size(EdB, 2)]) + 20 * log10(R ./ R0) + 2 * a * (R - R0);
    % 
    % Store the computed sy back into the structure or handle it as needed
    allData.(fields{i}).echo_intensity = E;

    % Optional: Display or store the results
    fprintf('Processed %s\n', fields{i});
end

% At this point, allData has a new field 'backscatter_strength' for each dataset
disp('All files have been processed and backscatter strength has been calculated.');
%% 

flds=fieldnames(allData);
ei=[];et=[];
figure; hold on
for i=1:length(flds)
    if allData.(flds{i}).instrument_depth>250
        
        [x,y]=meshgrid(allData.(flds{i}).distance,allData.(flds{i}).time);
        [xx,yy]=meshgrid(27,allData.(flds{i}).time);

        et=[et;allData.(flds{i}).time];
        ei=[ei; griddata(x,y,double(allData.(flds{i}).echo_intensity),xx,yy)];
       
    end
end
%%
figure
plot(et,smoothdata(ei,'movmean',24*7,'omitnan'));


load BAKmay2024.mat
yyaxis right
plot(BAK.mtime,smoothdata(BAK.B,'movmean',24*7,'omitnan'));
xlim([datenum(2017,01,01) datenum(2024,01,01)]);
axdate