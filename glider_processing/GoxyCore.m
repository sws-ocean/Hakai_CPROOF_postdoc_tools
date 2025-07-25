% Oxygen corrections

clear 
load glider_allProfs allProfs

ARODcal=rdAROD('AROD_cal.csv');
AADIcal=rdAADI('AADI_cal.csv');

allProfs.s_t=gsw_sigma0(allProfs.SA,allProfs.CT);

load DO_1990-2024_QCS_CTD.mat
QCSctd.s_t=gsw_sigma0(QCSctd.SA,QCSctd.CT);

load OS17-23_ctd.mat

%% Define mission list and run initial corrections

% find Calvert line missions
Mmsk=allProfs.missionIdx==1;
mission_list=unique(vertcat(allProfs.mission{Mmsk}));

mission_list(strcmp(mission_list,...
    'dfo-colin1142-20240212_grid_delayed.nc'))=[];

% Run through missions and apply T and P compensation

d = dir('/Users/samuelstevens/Dropbox/Hakai/gliders/data/delayed_ungrid/outdir/*nc');
tau=NaN(length(mission_list));
r=NaN(length(mission_list));
tic

allProfs.oxygen_corrected=NaN(size(allProfs.oxygen_concentration));

for i=1:length(mission_list)

    % index the correct mission
    M=mission_list{i};
    parts = strsplit(M, '-');
    fmsk=contains(string({d(:).name}),[parts{2},'-',parts{3}(1:8)]);

       % find mission in dataset and extract 
    Mmsk=contains(string(allProfs.mission),mission_list{i});
    op=unique([allProfs.optode{Mmsk}]);

    % Use regex to extract the serial number if the pattern is found
    serial_pattern = '''serial'':\s*''(\d+)''';
    serial_number = regexp(op, serial_pattern, 'tokens', 'once');

    if ~isempty(serial_number)
        serial_number = serial_number{1}; % Extract the serial number from the cell
    else
        disp('Serial number not found');
    end


    if contains(op,'AROD_FT') && sum(fmsk)>0
        finfo=ncinfo(d(fmsk).name);
        if sum(contains(string({finfo.Variables.Name}),'oxygen'))>0

            % Read in raw data
            rawO = ncread(d(fmsk).name, 'oxygen_concentration');
            % skip iteration if no O2 data available
            if sum(~isnan(rawO))==0
                continue;
            end
            rawS = ncread(d(fmsk).name, 'salinity'); % PSU
            rawT = ncread(d(fmsk).name, 'temperature');
            rawP = ncread(d(fmsk).name, 'pressure');
            rawL = ncread(d(fmsk).name, 'longitude');
            rawLL = ncread(d(fmsk).name, 'latitude');
            Pidx = ncread(d(fmsk).name, 'profile_index');
            rawTime = ncread(d(fmsk).name, 'time');
            rawT(rawT>40)=NaN;

            Bcoeffs={-0.00624523;-0.00737614;-0.01034100;-0.00817083;-0.000000488682}; % B0 B1 B2 B3 C0
            sal_ref = 0; % These are always the case for Rinko sensors

            % Find in cal sheet
            calMsk=str2double(serial_number)==ARODcal.SerialNumber;
            Cp=ARODcal.Cp(calMsk);

            % pressure correction
            p_factor = (rawP / 100) .* Cp + 1;

            % temperature correction
            Ts = log((298.15 - rawT) ./ (273.15 + rawT)); % scaled temperature
            s_factor = exp((rawS - sal_ref) .* (Bcoeffs{1} + Bcoeffs{2} .* Ts + Bcoeffs{3} .* Ts.^2 + Bcoeffs{4} .* Ts.^3) + Bcoeffs{5} .* (rawS.^2-sal_ref^2));

            % apply corrections to the uncorrected oxygen concentration
            corrO = rawO .* p_factor .* s_factor;

            % Phase calculation (Uchida et al., 2010)
            d0 = ARODcal.d0(calMsk); % from calibration data
            d1 = ARODcal.d1(calMsk);
            d2 = ARODcal.d2(calMsk);

            c0=ARODcal.c0(calMsk);
            c1=ARODcal.c1(calMsk);
            c2=ARODcal.c2(calMsk);

            % Calculate Ksv
            Ksv = c0 + (c1 .* rawT) + (c2 .* rawT.^2);

            % Calculate phase from corrected oxygen concentration
            oxygen_phase = (((1 + d0 .* rawT) ./ ((corrO .* Ksv) + 1)) - d1) ./ d2;

            pgrid = [0:1100]'; % Regular pressure grid

            uPhase = []; dcPhase = [];
            upDt = []; dnDt = [];
            upT=[]; dnT=[];

            for ii = 0:2:max(Pidx)
                dmsk = Pidx == ii & rawP > 10 & ~isnan(rawO);
                umsk = Pidx == ii+1 & rawP > 10 & ~isnan(rawO);

                if sum(dmsk) > 1 && sum(umsk) > 1

                    % interpolate onto regular depth grids
                    dnDt =[dnDt;interp1(rawP(dmsk)+rand(size(rawP(dmsk)))/1e6,...
                        rawTime(dmsk),pgrid)];

                    dcPhase =[dcPhase;interp1(rawP(dmsk)+rand(size(rawP(dmsk))),...
                        oxygen_phase(dmsk),pgrid)];

                    upDt =[upDt;interp1(rawP(umsk)+rand(size(rawP(umsk)))/1e6,...
                        rawTime(umsk),pgrid)];

                    uPhase =[uPhase;interp1(rawP(umsk)+rand(size(rawP(umsk))),...
                        oxygen_phase(umsk),pgrid)];

                    upT =[upT;interp1(rawP(umsk)+rand(size(rawP(umsk))),...
                        rawT(umsk),pgrid)];

                    dnT =[dnT;interp1(rawP(dmsk)+rand(size(rawP(dmsk))),...
                        rawT(dmsk),pgrid)];

                    % else
                    %     dnDt =[dnDt;NaN(length(pgrid),1)];
                    %     dcPhase =[dcPhase;NaN(length(pgrid),1)];
                    %     uPhase =[uPhase;NaN(length(pgrid),1)];
                    %     upDt =[upDt;NaN(length(pgrid),1)];
                    %     upT =[upT;NaN(length(pgrid),1)];
                    %     dnT =[dnT;NaN(length(pgrid),1)];

                end
            end

            % Remove outliers in dt
            up_dt = abs(gradient(upDt));
            down_dt = abs(gradient(dnDt));

            valid_up = up_dt < 30 & up_dt > 0.1;
            valid_down = down_dt < 30 & down_dt > 0.1;
            up_field = uPhase; up_field(~valid_up) = NaN;
            dn_field = dcPhase; dn_field(~valid_down) = NaN;

            % Calculate derivatives using finite difference with respect to time intervals
            valid_points = ~isnan(up_field) & ~isnan(dn_field);

            % Filter the fields to only valid points
            filtUpPhase = up_field(valid_points);
            filtDnPhase = dn_field(valid_points);
            filtUpT = upDt(valid_points);
            filtDnT = dnDt(valid_points);

            % Calculate centered finite difference for upcast and downcast phases
            dFdt_upcast = diff(filtUpPhase) ./ diff(filtUpT);
            dFdt_downcast = diff(filtDnPhase) ./ diff(filtDnT);

            % Now we calculate y and x for the regression
            % Exclude the last point from filtered fields to match the size of the derivative arrays
            y = filtUpPhase(1:end-1) - filtDnPhase(1:end-1);
            x = dFdt_upcast - dFdt_downcast;

            % Perform linear regression to find slope, which is -tau
            slope = polyfit(x, y, 1);
            tau(i) = -slope(1);

            % Print tau result
            tmp = corrcoef(x, y);
            r(i) = tmp(2);
            fprintf('Tau for AROD_FT %s is: %.4f seconds (r = %.4f)\n', serial_number, tau(i), r(i));
            
            % interpolate and boxcar filter phase %
            timeGrid=[rawTime(1):rawTime(end)]';
            vMsk = ~isnan(rawTime) & ~isnan(oxygen_phase);
            Phasei=interp1(rawTime(vMsk),oxygen_phase(vMsk),timeGrid);
            vMsk = ~isnan(rawTime) & ~isnan(rawP);
            Pi=interp1(rawTime(vMsk),rawP(vMsk),timeGrid);
            vMsk = ~isnan(rawTime) & ~isnan(Pidx);
            Pidxi=interp1(rawTime(vMsk),Pidx(vMsk),timeGrid,'nearest');
            vMsk = ~isnan(rawTime) & ~isnan(rawT);
            Ti=interp1(rawTime(vMsk),rawT(vMsk),timeGrid);

            smthPhase = movmean(Phasei, tau(i), 'omitnan', 'SamplePoints', timeGrid);
            dt=diff(timeGrid);

            % Initialize array for corrected phase
            phaseCorrected = NaN(size(smthPhase));

            % Calculate 'a' and 'b' for each time interval, then apply correction
            for j = 1:length(smthPhase)-1
                % Calculate b based on the time difference and tau
                b = (1 + 2 * (tau(i) / dt(j)))^(-1);
                a = 1 - 2 * b;

                % Apply the inverse bilinear transform for phase correction
                phaseCorrected(j) = (1 / (2 * b)) * (smthPhase(j+1) - a * smthPhase(j));
            end
            
            % Calculate Ksv
            Ksv = c0 + (c1 .* Ti) + (c2 .* Ti.^2);

            % Reverse phase calculation to get oxygen concentration from corrected phase
            oxygen_concentration_corrected = ((1 + d0 .* Ti) ./...
                (phaseCorrected .* d2 + d1) - 1) ./ Ksv;

            figure('color','w');
            subplot(1,2,1)
            x=allProfs.oxygen_concentration(:,Mmsk);
            y=allProfs.pressure(:,Mmsk);
            scatter(x(:),y(:),2,'k','filled');
            hold on
            scatter(corrO,rawP,2,'r','filled');
            axis ij; grid on;ylim([0 1100]);

            subplot(1,2,2)
            x=allProfs.oxygen_concentration(:,Mmsk);
            y=allProfs.pressure(:,Mmsk);
            s1=scatter(x(:),y(:),2,'k','filled');        hold on
            s2= scatter(corrO,rawP,2,'r','filled');
            s3=scatter(oxygen_concentration_corrected,Pi,2,'g','filled');
            axis ij; grid on;ylim([0 1100]);
            sgtitle([parts{2},'-',parts{3}(1:8)]);
            text(100,650,sprintf('AROD-FT (Serial=%s)\n tau = %.1f s\n r = %.4f',...
                serial_number, tau(i), r(i)),'interpreter','Latex');
            legend([s1 s2 s3],'raw','P & S compensated','Tau corrected','location','southeast');

            set(findall(gcf,'-property','fontsize'),'fontsize',10);
            set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
            drawnow

            % Define bin edges
            binEdges = 0.5:1:1100.5;
            idxs=find(Mmsk);
            % Loop through each profile
            for j = 1:length(idxs)
                % Mask for current profile
                profileMask = (Pidxi == j);

                % Get Pi and oxygen values for the current profile
                Pi_profile = Pi(profileMask);
                oxygen_profile = oxygen_concentration_corrected(profileMask);

                % Bin the Pi values to get the bin indices
                [~, ~, binIdx] = histcounts(Pi_profile, binEdges);

                % Calculate the mean oxygen concentration for each bin using accumarray
                oxygen_means = accumarray(binIdx(binIdx > 0),...
                    oxygen_profile(binIdx > 0), [1100, 1], @mean, NaN);

                % Store the results
                allProfs.oxygen_corrected(:, idxs(j)) = oxygen_means;
            end


            eval(['export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/code/quickCorrFigs/',[parts{2},'-',parts{3}(1:8)],'.png -m5 -nofontswap']);
            close
        end

    elseif contains(op,'AADI') && sum(fmsk)>0
        finfo=ncinfo(d(fmsk).name);
        if sum(contains(string({finfo.Variables.Name}),'oxygen'))>0

            % Read in raw data
            rawO = ncread(d(fmsk).name, 'oxygen_concentration');
            % skip iteration if no O2 data available
            if sum(~isnan(rawO))==0
                continue;
            end
            rawS = ncread(d(fmsk).name, 'salinity'); % PSU
            rawT = ncread(d(fmsk).name, 'temperature');
            rawP = ncread(d(fmsk).name, 'pressure');
            rawL = ncread(d(fmsk).name, 'longitude');
            rawLL = ncread(d(fmsk).name, 'latitude');
            Pidx = ncread(d(fmsk).name, 'profile_index');
            rawTime = ncread(d(fmsk).name, 'time');

            % Remove crappy profiles
            msk=rawO>(mean(rawO,'omitnan')+3*std(rawO,1,'omitnan')) | ...
                rawO<(mean(rawO,'omitnan')-3*std(rawO,1,'omitnan'));
            if sum(msk)>0
                p=unique(Pidx(msk));
                for j=1:length(p)
                    rawO(Pidx==p(j))=NaN;
                end
            end

            % Do annoying corrections:
            if sum(contains(mission_list{i},'marvin1003-20230929'))>0
                rawO(Pidx==35)=NaN;
                rawO(Pidx==0)=NaN;
            elseif sum(contains(mission_list{i},'hal1002-20220804'))>0
                rawO(Pidx==960)=NaN; 
                rawO(Pidx==964)=NaN;
                rawO(Pidx==965)=NaN;
            end
            rawT(rawT>40)=NaN;

            Bcoeffs = {-0.00624097; -0.00693498; -0.00690358; -0.00429155; -3.1168e-7}; % B0 B1 B2 B3 C0
            sal_ref = 0; % [PSU]

            % Find in cal sheet
            calMsk=str2double(serial_number)==AADIcal.SerialNumber;
            Cp=0.032;

            % pressure correction
            p_factor = (rawP / 100) .* Cp + 1;

            % temperature correction
            Ts = log((298.15 - rawT) ./ (273.15 + rawT)); % scaled temperature
            s_factor = exp((rawS - sal_ref) .* (Bcoeffs{1} + Bcoeffs{2} .* Ts + Bcoeffs{3} .* Ts.^2 + Bcoeffs{4} .* Ts.^3) + Bcoeffs{5} .* (rawS.^2-sal_ref^2));

            % apply corrections to the uncorrected oxygen concentration
            corrO = rawO .* p_factor .* s_factor;

            % add coefficients and calcularte phase here... %

            % Calculate phase (Uchida et al., 2008 modified Stern-Volmer equation)
            c0 = AADIcal.SVUfoilCoef_0(calMsk);
            c1 = AADIcal.SVUfoilCoef_1(calMsk);
            c2 = AADIcal.SVUfoilCoef_2(calMsk);
            c3 = AADIcal.SVUfoilCoef_3(calMsk);
            c4 = AADIcal.SVUfoilCoef_4(calMsk);
            c5 = AADIcal.SVUfoilCoef_5(calMsk);
            c6 = AADIcal.SVUfoilCoef_6(calMsk);

            % Compute Ksv and P0 using temperature
            Ksv = c0 + c1 .* rawT + c2 .* rawT.^2;
            P0 = c3 + c4 .* rawT;

            % % Calculate phase from compensated oxygen concentration
            oxygen_phase = ((P0 ./ ((corrO .* Ksv) + 1)) - c5) / c6;

            pgrid = [0:1100]'; % Regular pressure grid

            uPhase = []; dcPhase = [];
            upDt = []; dnDt = [];
            upT=[]; dnT=[];

            for ii = 0:2:max(Pidx)
                dmsk = Pidx == ii & rawP > 10 & ~isnan(rawO);
                umsk = Pidx == ii+1 & rawP > 10 & ~isnan(rawO);

                if sum(dmsk) > 1 && sum(umsk) > 1

                    % interpolate onto regular depth grids
                    dnDt =[dnDt;interp1(rawP(dmsk)+rand(size(rawP(dmsk)))/1e6,...
                        rawTime(dmsk),pgrid)];

                    dcPhase =[dcPhase;interp1(rawP(dmsk)+rand(size(rawP(dmsk))),...
                        oxygen_phase(dmsk),pgrid)];

                    upDt =[upDt;interp1(rawP(umsk)+rand(size(rawP(umsk)))/1e6,...
                        rawTime(umsk),pgrid)];

                    uPhase =[uPhase;interp1(rawP(umsk)+rand(size(rawP(umsk))),...
                        oxygen_phase(umsk),pgrid)];

                    upT =[upT;interp1(rawP(umsk)+rand(size(rawP(umsk))),...
                        rawT(umsk),pgrid)];

                    dnT =[dnT;interp1(rawP(dmsk)+rand(size(rawP(dmsk))),...
                        rawT(dmsk),pgrid)];

                    % else
                    %     dnDt =[dnDt;NaN(length(pgrid),1)];
                    %     dcPhase =[dcPhase;NaN(length(pgrid),1)];
                    %     uPhase =[uPhase;NaN(length(pgrid),1)];
                    %     upDt =[upDt;NaN(length(pgrid),1)];
                    %     upT =[upT;NaN(length(pgrid),1)];
                    %     dnT =[dnT;NaN(length(pgrid),1)];

                end
            end

            % Remove outliers in dt
            up_dt = abs(gradient(upDt));
            down_dt = abs(gradient(dnDt));

            valid_up = up_dt < 60 & up_dt > 0.1;
            valid_down = down_dt < 60 & down_dt > 0.1;
            up_field = uPhase; up_field(~valid_up) = NaN;
            dn_field = dcPhase; dn_field(~valid_down) = NaN;

            % Calculate derivatives using finite difference with respect to time intervals
            valid_points = ~isnan(up_field) & ~isnan(dn_field);

            % Filter the fields to only valid points
            filtUpPhase = up_field(valid_points);
            filtDnPhase = dn_field(valid_points);
            filtUpT = upDt(valid_points);
            filtDnT = dnDt(valid_points);

            % Calculate centered finite difference for upcast and downcast phases
            dFdt_upcast = diff(filtUpPhase) ./ diff(filtUpT);
            dFdt_downcast = diff(filtDnPhase) ./ diff(filtDnT);

            % Now we calculate y and x for the regression
            % Exclude the last point from filtered fields to match the size of the derivative arrays
            y = filtUpPhase(1:end-1) - filtDnPhase(1:end-1);
            x = dFdt_upcast - dFdt_downcast;

            % Perform linear regression to find slope, which is -tau
            slope = polyfit(x, y, 1);
            tau(i) = -slope(1);

            % Print tau result
            tmp = corrcoef(x, y);
            r(i) = tmp(2);
            fprintf('Tau for AADI %s is: %.4f seconds (r = %.4f)\n', serial_number, tau(i), r(i));

            % interpolate and boxcar filter phase %
            timeGrid=[rawTime(1):rawTime(end)]';
            vMsk = ~isnan(rawTime) & ~isnan(oxygen_phase);
            Phasei=interp1(rawTime(vMsk),oxygen_phase(vMsk),timeGrid);
            vMsk = ~isnan(rawTime) & ~isnan(rawP);
            Pi=interp1(rawTime(vMsk),rawP(vMsk),timeGrid);
            vMsk = ~isnan(rawTime) & ~isnan(Pidx);
            Pidxi=interp1(rawTime(vMsk),Pidx(vMsk),timeGrid,'nearest');
            vMsk = ~isnan(rawTime) & ~isnan(rawT);
            Ti=interp1(rawTime(vMsk),rawT(vMsk),timeGrid);

            smthPhase = movmean(Phasei, tau(i), 'omitnan', 'SamplePoints', timeGrid);
            dt=diff(timeGrid);

            % Initialize array for corrected phase
            phaseCorrected = NaN(size(smthPhase));

            % Calculate 'a' and 'b' for each time interval, then apply correction
            for j = 1:length(smthPhase)-1
                % Calculate b based on the time difference and tau
                b = (1 + 2 * (tau(i) / dt(j)))^(-1);
                a = 1 - 2 * b;

                % Apply the inverse bilinear transform for phase correction
                phaseCorrected(j) = (1 / (2 * b)) * (smthPhase(j+1) - a * smthPhase(j));
            end

             % Compute Ksv and P0 using temperature
            Ksv = c0 + c1 .* Ti + c2 .* Ti.^2;
            P0 = c3 + c4 .* Ti;

             % Reverse the calculation to get oxygen concentration from corrected phase
            oxygen_concentration_corrected = (P0 ./...
                ((phaseCorrected .* c6 + c5)) - 1) ./ Ksv;

            figure('color','w');
            subplot(1,2,1)
            x=allProfs.oxygen_concentration(:,Mmsk);
            y=allProfs.pressure(:,Mmsk);
            scatter(x(:),y(:),2,'k','filled');
            hold on
            scatter(corrO,rawP,2,'r','filled');
            axis ij; grid on;ylim([0 1100]);

            subplot(1,2,2)
            x=allProfs.oxygen_concentration(:,Mmsk);
            y=allProfs.pressure(:,Mmsk);
            s1=scatter(x(:),y(:),2,'k','filled');        hold on
            s2= scatter(corrO,rawP,2,'r','filled');
            s3=scatter(oxygen_concentration_corrected,Pi,2,'g','filled');
            axis ij; grid on;ylim([0 1100]);
            sgtitle([parts{2},'-',parts{3}(1:8)]);
            text(100,650,sprintf('AADI (Serial=%s)\n tau = %.1f s\n r = %.4f',...
                serial_number, tau(i), r(i)),'interpreter','Latex');
            legend([s1 s2 s3],'raw','P & S compensated','Tau corrected','location','southeast');

            set(findall(gcf,'-property','fontsize'),'fontsize',10);
            set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
            drawnow

            % Define bin edges
            binEdges = 0.5:1:1100.5;
            idxs=find(Mmsk);

            % Loop through each profile
            for j = 1:length(idxs)
                % Mask for current profile
                profileMask = (Pidxi == j);

                % Get Pi and oxygen values for the current profile
                Pi_profile = Pi(profileMask);
                oxygen_profile = oxygen_concentration_corrected(profileMask);

                % Bin the Pi values to get the bin indices
                [~, ~, binIdx] = histcounts(Pi_profile, binEdges);

                % Calculate the mean oxygen concentration for each bin using accumarray
                oxygen_means = accumarray(binIdx(binIdx > 0),...
                    oxygen_profile(binIdx > 0), [1100, 1], @mean, NaN);

                % Store the results
                allProfs.oxygen_corrected(:, idxs(j)) = oxygen_means;
            end

            eval(['export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/code/quickCorrFigs/',[parts{2},'-',parts{3}(1:8)],'.png -m5 -nofontswap']);
            close
        end

    elseif contains(op,'JFE') && sum(fmsk)>0 && sum(contains(string({finfo.Variables.Name}),'oxygen'))>0
        finfo=ncinfo(d(fmsk).name);
        if sum(contains(string({finfo.Variables.Name}),'oxygen'))>0

            % Read in raw data
            rawOs = ncread(d(fmsk).name, 'oxygen_saturation');
            rawS = ncread(d(fmsk).name, 'salinity'); % PSU
            rawT = ncread(d(fmsk).name, 'temperature');
            rawP = ncread(d(fmsk).name, 'pressure');
            rawL = ncread(d(fmsk).name, 'longitude');
            rawLL = ncread(d(fmsk).name, 'latitude');
            Pidx = ncread(d(fmsk).name, 'profile_index');
            rawPD=ncread(d(fmsk).name, 'potential_density');
            rawTime = ncread(d(fmsk).name, 'time');

            % Only one AROD-CAR sensors each with a single coeff
            E=0.00440000;

            % pressure correction
            corrO= rawOs.*(1 + E.*rawP.*0.01);

            % convert saturation % to concentration
            SA=gsw_SA_from_SP(rawS,rawP,rawL,rawLL);
            CT=gsw_CT_from_t(SA,rawT,rawP);

            OSat = gsw_O2sol(SA,CT,rawP,rawL,rawLL);
            corrO = OSat.*corrO.*rawPD./(100*1000);

            pgrid = [0:1100]'; % Regular pressure grid

            uPhase = []; dcPhase = [];
            upDt = []; dnDt = [];

            for ii = 0:2:max(Pidx)
                dmsk = Pidx == ii & rawP > 10 & ~isnan(rawOs);
                umsk = Pidx == ii+1 & rawP > 10 & ~isnan(rawOs);

                if sum(dmsk) > 1 && sum(umsk) > 1

                    % interpolate onto regular depth grids
                    dnDt =[dnDt;interp1(rawP(dmsk)+rand(size(rawP(dmsk)))/1e6,...
                        rawTime(dmsk),pgrid)];

                    dcPhase =[dcPhase;interp1(rawP(dmsk)+rand(size(rawP(dmsk))),...
                        corrO(dmsk),pgrid)];

                    upDt =[upDt;interp1(rawP(umsk)+rand(size(rawP(umsk)))/1e6,...
                        rawTime(umsk),pgrid)];

                    uPhase =[uPhase;interp1(rawP(umsk)+rand(size(rawP(umsk))),...
                        corrO(umsk),pgrid)];

                end
            end

            % Remove outliers in dt
            up_dt = abs(gradient(upDt));
            down_dt = abs(gradient(dnDt));

            valid_up = up_dt < 30 & up_dt > 0.1;
            valid_down = down_dt < 30 & down_dt > 0.1;
            up_field = uPhase; up_field(~valid_up) = NaN;
            dn_field = dcPhase; dn_field(~valid_down) = NaN;

            % Calculate derivatives using finite difference with respect to time intervals
            valid_points = ~isnan(up_field) & ~isnan(dn_field);

            % Filter the fields to only valid points
            filtUpPhase = up_field(valid_points);
            filtDnPhase = dn_field(valid_points);
            filtUpT = upDt(valid_points);
            filtDnT = dnDt(valid_points);

            % Calculate centered finite difference for upcast and downcast phases
            dFdt_upcast = diff(filtUpPhase) ./ diff(filtUpT);
            dFdt_downcast = diff(filtDnPhase) ./ diff(filtDnT);

            % Now we calculate y and x for the regression
            % Exclude the last point from filtered fields to match the size of the derivative arrays
            y = filtUpPhase(1:end-1) - filtDnPhase(1:end-1);
            x = dFdt_upcast - dFdt_downcast;

            % Perform linear regression to find slope, which is -tau
            slope = polyfit(x, y, 1);
            tau(i) = -slope(1);

            % Print tau result
            tmp = corrcoef(x, y);
            r(i) = tmp(2);
            fprintf('Tau for ARO-CAR-Z10 %s is: %.4f seconds (r = %.4f)\n', serial_number, tau(i), r(i));

            % interpolate and boxcar filter phase %
            timeGrid=[rawTime(1):rawTime(end)]';
            vMsk = ~isnan(rawTime) & ~isnan(corrO);
            Phasei=interp1(rawTime(vMsk),corrO(vMsk),timeGrid);
            vMsk = ~isnan(rawTime) & ~isnan(rawP);
            Pi=interp1(rawTime(vMsk),rawP(vMsk),timeGrid);
            vMsk = ~isnan(rawTime) & ~isnan(Pidx);
            Pidxi=interp1(rawTime(vMsk),Pidx(vMsk),timeGrid,'nearest');
            vMsk = ~isnan(rawTime) & ~isnan(rawT);
            Ti=interp1(rawTime(vMsk),rawT(vMsk),timeGrid);

            smthPhase = movmean(Phasei, tau(i), 'omitnan', 'SamplePoints', timeGrid);
            dt=diff(timeGrid);

            % Initialize array for corrected phase
            oxygen_concentration_corrected = NaN(size(smthPhase));

            % Calculate 'a' and 'b' for each time interval, then apply correction
            for j = 1:length(smthPhase)-1
                % Calculate b based on the time difference and tau
                b = (1 + 2 * (tau(i) / dt(j)))^(-1);
                a = 1 - 2 * b;

                % Apply the inverse bilinear transform for phase correction
                oxygen_concentration_corrected(j) = (1 / (2 * b)) * (smthPhase(j+1) - a * smthPhase(j));
            end

            figure('color','w');
            subplot(1,2,1)
            x=allProfs.oxygen_concentration(:,Mmsk);
            y=allProfs.pressure(:,Mmsk);
            scatter(x(:),y(:),2,'k','filled');
            hold on
            scatter(corrO,rawP,2,'r','filled');
            axis ij; grid on;ylim([0 1100]);

            subplot(1,2,2)
            x=allProfs.oxygen_concentration(:,Mmsk);
            y=allProfs.pressure(:,Mmsk);
            s1=scatter(x(:),y(:),2,'k','filled');        hold on
            s2= scatter(corrO,rawP,2,'r','filled');
            s3=scatter(oxygen_concentration_corrected,Pi,2,'g','filled');
            axis ij; grid on;ylim([0 1100]);
            sgtitle([parts{2},'-',parts{3}(1:8)]);
            text(100,650,sprintf('ARO-CAR-Z10 (Serial=%s)\n tau = %.1f s\n r = %.4f',...
                serial_number, tau(i), r(i)),'interpreter','Latex');
            legend([s1 s2 s3],'raw','P compensated','Tau corrected','location','southeast');

            set(findall(gcf,'-property','fontsize'),'fontsize',10);
            set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
            drawnow

            eval(['export_fig /Users/samuelstevens/Dropbox/Hakai/gliders/oxygen/code/quickCorrFigs/',[parts{2},'-',parts{3}(1:8)],'.png -m5 -nofontswap']);
            close

            % Define bin edges
            binEdges = 0.5:1:1100.5;
            idxs=find(Mmsk);

            % Loop through each profile
            for j = 1:length(idxs)
                % Mask for current profile
                profileMask = (Pidxi == j);

                % Get Pi and oxygen values for the current profile
                Pi_profile = Pi(profileMask);
                oxygen_profile = oxygen_concentration_corrected(profileMask);

                % Bin the Pi values to get the bin indices
                [~, ~, binIdx] = histcounts(Pi_profile, binEdges);

                % Calculate the mean oxygen concentration for each bin using accumarray
                oxygen_means = accumarray(binIdx(binIdx > 0),...
                    oxygen_profile(binIdx > 0), [1100, 1], @mean, NaN);

                % Store the results
                allProfs.oxygen_corrected(:, idxs(j)) = oxygen_means;
            end

        end
    end
end

toc

save('/Users/samuelstevens/Dropbox/Hakai/gliders/data/corrected_allProfs_May25.mat','allProfs',"-v7.3");
close all;


%% Compare with CTD and apply drift corrections
clear
load corrected_allProfs.mat
load CTDcompv2.mat

mission_list={'dfo-eva035-20230518','dfo-eva035-20230620',...
    'dfo-eva035-20230720','dfo-eva035-20230811','dfo-eva035-20230915'...
    'dfo-bb046-20220507','dfo-bb046-20220608','dfo-hal1002-20220804',...
    'dfo-marvin1003-20221018','dfo-marvin1003-20220928'};

flds=fieldnames(CTDcomp);

sigGrid=[20:0.05:27]';
sigIdx=find(sigGrid==25);

% convert CTD ml/l to umol/kg and calculate sigma_theta
for i=1:length(flds)
    for ii=1:length(CTDcomp.(flds{i}).hakai_id)
        CTDcomp.(flds{i}).SA{ii}=gsw_SA_from_SP(CTDcomp.(flds{i}).salinity{ii},...
            double(CTDcomp.(flds{i}).pressure{ii}),CTDcomp.(flds{i}).longitude{ii},...
            CTDcomp.(flds{i}).latitude{ii});
        
        CTDcomp.(flds{i}).CT{ii}=gsw_CT_from_t(CTDcomp.(flds{i}).SA{ii},...
            CTDcomp.(flds{i}).temperature{ii},double(CTDcomp.(flds{i}).pressure{ii}));

        CTDcomp.(flds{i}).rho{ii}=gsw_rho(CTDcomp.(flds{i}).SA{ii},...
            CTDcomp.(flds{i}).CT{ii},double(CTDcomp.(flds{i}).pressure{ii}));

        CTDcomp.(flds{i}).dissolved_oxygen_ml_l{ii}=...
            (CTDcomp.(flds{i}).dissolved_oxygen_ml_l{ii} * 44660) ./ ...
            CTDcomp.(flds{i}).rho{ii};

        CTDcomp.(flds{i}).s_t{ii}=gsw_sigma0(CTDcomp.(flds{i}).SA{ii},...
            CTDcomp.(flds{i}).CT{ii});
    end
end

err=[];
for i=1:length(mission_list)
    Mmsk=contains(string(allProfs.mission),mission_list{i});

    uncorrO=allProfs.oxygen_concentration(:,Mmsk);
    pres=allProfs.pressure(:,Mmsk);
    sal=allProfs.salinity(:,Mmsk);
    temp=allProfs.temperature(:,Mmsk);
    lat=allProfs.latitude(Mmsk);
    lon=allProfs.longitude(Mmsk);
    time=allProfs.mtime(Mmsk);
    s_t=allProfs.s_t(:,Mmsk);
    corrO = allProfs.oxygen_corrected(:,Mmsk);

    % Compare with CTD profiles
    % index the correct mission
    M=mission_list{i};
    parts = strsplit(M, '-');
    msk=contains(string(flds),parts{3});

    if sum(msk)>0
        figure('color','w');
        coeffs=[];

        for ii=1:length(CTDcomp.(flds{msk}).hakai_id)

            subplot(4,length(CTDcomp.(flds{msk}).hakai_id),ii);
            hold on;

            ddx = gsw_distance([lon' ones(size(lon))' .* CTDcomp.(flds{msk}).longitude{ii}(1)], ...
                [lat' ones(size(lat))' .* CTDcomp.(flds{msk}).latitude{ii}(1)]);
            ddt = abs(CTDcomp.(flds{msk}).time{ii}(1) - time);
            % [dx(ii),~]=min(ddx);

            msk1=find(ddx'<10000 & max(pres)>125 & sum(~isnan(s_t))>0 &...
                sum(~isnan(corrO))>0);

            if ~isempty(msk1)

                [dt,idx]=min(ddt(msk1));

                % interpolate onto regular density grid
                v=~isnan(s_t(:,msk1(idx)));
                Gco=interp1(s_t(v,msk1(idx)),corrO(v,msk1(idx)),sigGrid);
                Guo=interp1(s_t(v,msk1(idx)),uncorrO(v,msk1(idx)),sigGrid);
                CTDo=interp1(CTDcomp.(flds{msk}).s_t{ii},...
                    CTDcomp.(flds{msk}).dissolved_oxygen_ml_l{ii},sigGrid);

                plot(uncorrO(:,msk1(idx)),s_t(:,msk1(idx)),'r');
                plot(CTDcomp.(flds{msk}).dissolved_oxygen_ml_l{ii},...
                    CTDcomp.(flds{msk}).s_t{ii},'k');
                axis ij;grid on; xlim([25 350]);
                ddc=mean(Guo(sigIdx:end)-CTDo(sigIdx:end),'omitnan');
                sdc=mean(Guo(1:sigIdx)-CTDo(1:sigIdx),'omitnan');
                text(0.4,0.2,sprintf('MAD>%.2f kg/m3\n%3.0f umol/kg\n MAD<%.2f kg/m3\n%3.0f umol/kg\n  dx=%3.0f m\ndt=%3.0f mins',...
                    sigGrid(sigIdx),ddc,sigGrid(sigIdx),sdc, ddx(msk1(idx)),dt*24*60),'fontsize',7,'units','normalized');

                if ii==length(CTDcomp.(flds{msk}).hakai_id)
                    legend('Glider','CTD');
                end

                if ii==2
                    title('Uncorrected glider and CTD');
                end

                xlabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
                ylabel('\sigma_\theta (kg/m^3)');

                subplot(4,length(CTDcomp.(flds{msk}).hakai_id),ii+length(CTDcomp.(flds{msk}).hakai_id));hold on;
                plot(corrO(:,msk1(idx)),s_t(:,msk1(idx)),'r');
                plot(CTDcomp.(flds{msk}).dissolved_oxygen_ml_l{ii},...
                    CTDcomp.(flds{msk}).s_t{ii},'k');
                axis ij;grid on;
                ddc=mean(abs(Gco(sigIdx:end)-CTDo(sigIdx:end)),'omitnan');
                sc=mean(abs(Gco(1:sigIdx)-CTDo(1:sigIdx)),'omitnan');
                text(0.4,0.2,sprintf('MAD>%.2f kg/m3\n%3.1f umol/kg\n MAD<%.2f kg/m3\n%3.1f umol/kg\n',...
                    sigGrid(sigIdx),ddc,sigGrid(sigIdx),sc),'fontsize',7,'units','normalized');
                xlim([25 350]);
                xlabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
                ylabel('\sigma_\theta (kg/m^3)');

                if ii==2
                    title('P and S compensated glider and CTD');
                end

                % look at "drift"
                % Clean data using rmmissing
                clean_O= rmmissing([Gco(sigIdx:end),CTDo(sigIdx:end)]);
                try
                    % Linear least-squares fit (ignores NaNs automatically)
                    coeffs(ii) = clean_O(:,1) \ clean_O(:,2); % This finds the slope directly
                    subplot(4,length(CTDcomp.(flds{msk}).hakai_id),ii+2*length(CTDcomp.(flds{msk}).hakai_id));
                    hold on;
                    scatter(clean_O(:,1),clean_O(:,2) ,10);
                    ylabel('CTD DO (umol/kg)');
                    grid on;
                    plot(min(clean_O(:)):10:max(clean_O(:)),min(clean_O(:)):10:max(clean_O(:)),...
                        'color','k');

                    xlim([min(clean_O(:)) max(clean_O(:))]);
                    ylim([min(clean_O(:)) max(clean_O(:))]);
                    xticks(round(min(clean_O(:))):20:max(clean_O(:)));
                    yticks(round(min(clean_O(:))):20:max(clean_O(:)));
                    xlabel('Glider DO (umol/kg)');
                    if ii==2
                        title('P and S compensated glider vs CTD');
                    end
                    text(0.6,0.2,sprintf('m=%2.2f',coeffs(ii)),...
                        'fontsize',7,'units','normalized');

                end

            end
        end

        sgtitle(M);
        coeffs(coeffs==0)=[];

        dcorrO=mean(coeffs).*corrO;
        tmpO=mean(coeffs).*Gco;

        fprintf('Slope=%2.2f\n',max(coeffs));

        for ii=1:length(CTDcomp.(flds{msk}).hakai_id)

            ddx = gsw_distance([lon' ones(size(lon))' .* CTDcomp.(flds{msk}).longitude{ii}(1)], ...
                [lat' ones(size(lat))' .* CTDcomp.(flds{msk}).latitude{ii}(1)]);
            ddt = abs(CTDcomp.(flds{msk}).time{ii}(1) - time);
            % [dx(ii),~]=min(ddx);

            msk1=find(ddx'<10000 & max(pres)>125 & sum(~isnan(s_t))>0 &...
                sum(~isnan(corrO))>0);

            if ~isempty(msk1)
                [dt,idx]=min(ddt(msk1));

                % interpolate onto regular density grid
                v=~isnan(s_t(:,msk1(idx)));
                Gco=interp1(s_t(v,msk1(idx)),corrO(v,msk1(idx)),sigGrid);
                Guo=interp1(s_t(v,msk1(idx)),uncorrO(v,msk1(idx)),sigGrid);
                CTDo=interp1(CTDcomp.(flds{msk}).s_t{ii},...
                    CTDcomp.(flds{msk}).dissolved_oxygen_ml_l{ii},sigGrid);

                subplot(4,length(CTDcomp.(flds{msk}).hakai_id),ii+3*length(CTDcomp.(flds{msk}).hakai_id));hold on;
                ddx = gsw_distance([lon' ones(size(lon))' .* CTDcomp.(flds{msk}).longitude{ii}(1)], ...
                    [lat' ones(size(lat))' .* CTDcomp.(flds{msk}).latitude{ii}(1)]);
                ddt = abs(CTDcomp.(flds{msk}).time{ii}(1) - time);
                % [dx(ii),~]=min(ddx);

                msk1=find(ddx'<5000 & max(pres)>125);
                [dt,idx]=min(ddt(msk1));
                plot(dcorrO(:,msk1(idx)),s_t(:,msk1(idx)),'r');
                plot(CTDcomp.(flds{msk}).dissolved_oxygen_ml_l{ii},...
                    CTDcomp.(flds{msk}).s_t{ii},'k');
                axis ij;grid on;
                dc=mean(abs(tmpO(sigIdx:end)-CTDo(sigIdx:end)),'omitnan');
                sc=mean(abs(tmpO(1:sigIdx)-CTDo(1:sigIdx)),'omitnan');
                text(0.4,0.2,sprintf('MAD>%.2f kg/m3\n%3.1f umol/kg\n MAD<%.2f kg/m3\n%3.1f umol/kg\n',...
                    sigGrid(sigIdx),dc,sigGrid(sigIdx),sc),'fontsize',7,'units','normalized');
                xlim([25 350]);
                xlabel('$O$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');
                ylabel('\sigma_\theta (kg/m^3)');

                if ii==2
                    title(sprintf('Drift corrected (m=%2.2f) and CTD',mean(coeffs)));
                end

            end
        end

    else
        disp(['No comparison casts for ',M]);
        % fprintf('Using previous mission coeffs for this mision (slope=%2.2f)\n',max(coeffs));
        dcorrO=corrO;
    end

    a=input('Apply drift correction? y/n: ','s');
    if strcmp(a,'y')
        allProfs.oxygen_corrected(:,Mmsk)=dcorrO;
        err(i)=input('Error?: ');
    else
        err(i)=input('Error?: ');
    end

    set(findall(gcf,'-property','fontsize'),'fontsize',6);
    set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


    % % independent comparison with DFO data
    % for ii=1:length(QCSctd.latitude)
    %     ddx = gsw_distance([lon' ones(size(lon))' .* QCSctd.longitude(ii)], ...
    %         [lat' ones(size(lat))' .* QCSctd.latitude(ii)]);
    %     ddt = abs(QCSctd.time(ii) - time);
    %
    %     msk1=find(ddx'<10e3 & ddt < 10);
    %
    %     if ~isempty(msk1)
    %         % fprintf('x=%3.1f, t=%2.1f\n',min(ddx(msk1)),min(ddt(msk1)));
    %         [~,idx]=min(ddt(msk1));
    %
    %         figure;
    %         subplot(1,2,1)
    %         plot(dcorrO(:,msk1(idx)),pres(:,msk1(idx)),'r');
    %         hold on
    %         plot(corrO(:,msk1(idx)),pres(:,msk1(idx)),'g');
    %         plot(QCSctd.oxygen(:,ii),1:400,'k');
    %         grid on
    %         axis ij
    %
    %         text(0.6,0.2,sprintf('dx=%3.0f m\ndt=%3.1f days',ddx(msk1(idx)),...
    %             ddt(msk1(idx))),'fontsize',10,'units','normalized');
    %
    %         subplot(1,2,2)
    %         plot(dcorrO(:,msk1(idx)),s_t(:,msk1(idx)),'r');
    %         hold on
    %         plot(corrO(:,msk1(idx)),s_t(:,msk1(idx)),'g');
    %         plot(QCSctd.oxygen(:,ii),QCSctd.s_t(:,ii),'k');
    %         grid on
    %         axis ij
    %
    %     end
    %
    % end
    %
    % for ii=1:length(OSctd.lat)
    %
    %     ddx = gsw_distance([lon' ones(size(lon))' .* OSctd.lon(ii)], ...
    %         [lat' ones(size(lat))' .* OSctd.lat(ii)]);
    %     ddt = abs(OSctd.mtime(ii) - time);
    %
    %     msk1=find(ddx'<5e3 & ddt < 5);
    %
    %     if ~isempty(msk1)
    %         [~,idx]=min(ddt(msk1));
    %         % fprintf('x=%3.1f, t=%2.1f\n',ddx(msk1(idx)),ddt(msk1(idx)));
    %
    %         figure;
    %         subplot(1,2,1)
    %         plot(dcorrO(:,msk1(idx)),pres(:,msk1(idx)),'r');
    %         hold on
    %         plot(dcorrO(:,msk1(idx)),pres(:,msk1(idx)),'g');
    %         plot(OSctd.ox(:,ii),OSctd.pr(:,ii),'k');
    %         grid on
    %         axis ij
    %
    %         text(0.6,0.2,sprintf('dx=%3.0f m\ndt=%3.1f days',ddx(msk1(idx)),...
    %             ddt(msk1(idx))),'fontsize',10,'units','normalized');
    %
    %         subplot(1,2,2)
    %         plot(corrO(:,msk1(idx)),s_t(:,msk1(idx)),'r');
    %         hold on
    %         plot(dcorrO(:,msk1(idx)),s_t(:,msk1(idx)),'g');
    %         plot(OSctd.ox(:,ii),OSctd.s_t(:,ii),'k');
    %         grid on
    %         axis ij
    %     end
    % end
end

% save('/Users/samuelstevens/Dropbox/Hakai/gliders/data/Dcorrected_allProfs.mat','allProfs','err',"-v7.3");
close all;

