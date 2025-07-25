function data = readMooringCTDData(filename)
    % readMooringCTDData reads specified variables from a NetCDF file
    % associated with mooring CTD data, with flexibility for temperature and dissolved oxygen variable names.
    %
    % Usage:
    %   data = readMooringCTDData('path_to_file.nc')
    %
    % Input:
    %   filename - String or character array representing the path to the
    %              NetCDF file.
    %
    % Output:
    %   data - Struct containing the mooring data variables such as time,
    %          temperature, salinity, oxygen, pressure, latitude, and longitude.

    % Check if the NetCDF file exists
    if ~isfile(filename)
        error('File does not exist: %s', filename);
    end

    % Open the NetCDF file
    ncid = netcdf.open(filename, 'NC_NOWRITE');

    try
        % Initialize the data structure
        data = struct();
         
        % Read the time variable
        data.time = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'time'));

        % Read the temperature variable, checking for both potential names
        if hasVar(ncid, 'TEMPS901')
            data.temperature = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'TEMPS901'));
        elseif hasVar(ncid, 'TEMPST01')
            data.temperature = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'TEMPST01'));
        else
            error('Temperature variable not found under known names (TEMPS901 or TEMPST01)');
        end

        % Read salinity
        data.salinity = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'PSALST01'));

        % Read pressure
        data.pressure = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'PRESPR01'));
        
        % Check for oxygen data, or create a dummy NaN array if absent
        if hasVar(ncid, 'DOXYZZ01')
            data.oxygen = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'DOXYZZ01'));
        else
            % Create a NaN array of the same size as the temperature array
            data.oxygen = NaN(size(data.temperature), 'like', data.temperature);
        end

        % Read latitude and longitude
        if hasVar(ncid, 'latitude')
            data.latitude = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'latitude'));
        else
            error('Latitude variable not found.');
        end

        if hasVar(ncid, 'longitude')
            data.longitude = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'longitude'));
        else
            error('Longitude variable not found.');
        end

        % Calculate TEOS-10
        data.SA=gsw_SA_from_SP(data.salinity,data.pressure,data.longitude,...
            data.latitude);
        data.CT=gsw_CT_from_t(data.SA,data.temperature,data.pressure);
        data.sigma_theta=gsw_sigma0(data.SA,data.CT);
        data.rho=gsw_rho(data.SA,data.CT,data.pressure);
        data.spice=gsw_spiciness0(data.SA,data.CT);

        % convert ml/l to umol/kg
        data.oxygen=data.oxygen.* 44.661 ./ (data.rho/1000);

        % Calculate AOU
        data.PT=gsw_pt_from_CT(data.SA,data.CT);
        data.AOU=aou(data.salinity,data.PT,data.oxygen);

    catch ME
        netcdf.close(ncid);
        rethrow(ME);
    end

    % Close the NetCDF file
    netcdf.close(ncid);

    % Convert time to MATLAB datenum format
    data.time = datenum(1970, 1, 1, 0, 0, double(data.time));
end

function found = hasVar(ncid, varname)
    % Helper function to check if a variable exists in the netCDF file
    try
        netcdf.inqVarID(ncid, varname);
        found = true;
    catch
        found = false;
    end
end
