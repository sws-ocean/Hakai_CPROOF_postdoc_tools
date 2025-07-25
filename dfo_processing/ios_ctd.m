% SWS_ios_rd

% This script reads in all of the data from the IOS CTD archive,
% interpolates it onto a regular grid and stores it.
clear
cd /Users/samuelstevens/Dropbox/Hakai/ctd/data/hecate/outdir/
clc

%store all directories into a struct. .CHE= bottle files, .ctd=ctd files
ctd_dir.filnam1=dir('*.ctd'); ctd_dir.filnam1=dir('*.che'); 

%store all directory names in a struct
my_fields=fieldnames(ctd_dir); 

IOS_ctd.temp=NaN(length(1:1000),length(ctd_dir.filnam1));
IOS_ctd.sal=NaN(length(1:1000),length(ctd_dir.filnam1));
IOS_ctd.sal2=NaN(length(1:1000),length(ctd_dir.filnam1));
IOS_ctd.ox=NaN(length(1:1000),length(ctd_dir.filnam1));
% IOS_ctd.fluor=NaN(length(1:1000),length(ctd_dir.filnam1));
% IOS_ctd.no3=NaN(length(1:1000),length(ctd_dir.filnam1));
% IOS_ctd.chl=NaN(length(1:1000),length(ctd_dir.filnam1));
% IOS_ctd.xmiss=NaN(length(1:1000),length(ctd_dir.filnam1));
IOS_ctd.pr=repmat([1:1000]',1,length(ctd_dir.filnam1));
IOS_ctd.mtime=NaN(1,length(ctd_dir.filnam1));
IOS_ctd.lat=NaN(1,length(ctd_dir.filnam1));
IOS_ctd.lon=NaN(1,length(ctd_dir.filnam1));
IOS_ctd.station=cell(1,length(ctd_dir.filnam1));

ctd_count=0;

% Loop to extract and interpolate data
for i=1:length(my_fields) % loop through CTD directories
    
    for ii=1:length(ctd_dir.(my_fields{i})) % loop through all CTDs in directory
                
        filename=[ctd_dir.(my_fields{i})(ii).folder '/' ctd_dir.(my_fields{i})(ii).name];
        
        ctd=ios_rd(filename); % RPs IOS ctd read program
        
        if ~isfield(ctd,'pr')
            ctd.pr=gsw_p_from_z(ctd.depth*-1,ctd.latitude);
        end
       
        % If pressure nans are present, remove all data from that level.
        ctd.pr(isnan(ctd.pr))=[];
        
        % Depths need to be unique for interpolations
        [tmp,idx]=unique(ctd.pr);
        
        % Some casts are very shallow... need enough depth points for
        % interpolations. 
        if length(tmp)>2
            ctd_count=ctd_count+1;

            IOS_ctd.mtime(ctd_count)=ctd.mtime;
            IOS_ctd.lat(ctd_count)=ctd.latitude;
            IOS_ctd.lon(ctd_count)=ctd.longitude;
            
            if isfield(ctd,'temp')
                ctd.temp(isnan(ctd.pr))=[];
                IOS_ctd.temp(:,ctd_count)=interp1(ctd.pr(idx),ctd.temp(idx),[1:1000]');
            else
                ctd.temp2(isnan(ctd.pr))=[];
                IOS_ctd.temp(:,ctd_count)=interp1(ctd.pr(idx),ctd.temp2(idx),[1:1000]');
            end
            
            if isfield(ctd,'sal')
                ctd.sal(isnan(ctd.pr))=[];
                IOS_ctd.sal(:,ctd_count)=interp1(ctd.pr(idx),ctd.sal(idx),1:1000);
            end
            
            if isfield(ctd,'sal2')
                ctd.sal2(isnan(ctd.pr))=[];
                IOS_ctd.sal2(:,ctd_count)=interp1(ctd.pr(idx),ctd.sal2(idx),1:1000);
            end
            
            if isfield(ctd,'ox')
                % often there is ml/l and umol/l here. Only want to take ml/l.
                if isfield(ctd,'ox2')
                    if nanmean(ctd.ox) < nanmean(ctd.ox2)
                        ctd.ox(isnan(ctd.pr))=[];
                        IOS_ctd.ox(:,ctd_count)=interp1(ctd.pr(idx),ctd.ox(idx),1:1000);
                    else
                        ctd.ox2(isnan(ctd.pr))=[];
                        IOS_ctd.ox(:,ctd_count)=interp1(ctd.pr(idx),ctd.ox2(idx),1:1000);
                    end
                end
            end
            
            if isfield(ctd,'oxML')
                IOS_ctd.ox(:,ctd_count)=interp1(ctd.pr(idx),ctd.oxML(idx),1:1000);
            end
            
            
            if isfield(ctd,'station')
                IOS_ctd.station{ctd_count}=['IOS_' ctd.station];
            end
            
            
            
%             if isfield(ctd,'fluor')
%                 IOS_ctd.fluor(:,ctd_count)=interp1(ctd.pr(idx),ctd.fluor(idx),1:1000);
%             end
%             
%             if isfield(ctd,'no3')
%                 IOS_ctd.no3(:,ctd_count)=interp1(ctd.pr(idx),ctd.no3(idx),1:1000);
%             end
%             
%             if isfield(ctd,'chl')
%                 IOS_ctd.chl(:,ctd_count)=interp1(ctd.pr(idx),ctd.chl(idx),1:1000);
%             end
%             
%             if isfield(ctd,'xmiss')
%                 IOS_ctd.xmiss(:,ctd_count)=interp1(ctd.pr(idx),ctd.xmiss(idx),1:1000);
%             end
        end
    end
end

% Convert properties
IOS_ctd.sal=gsw_SR_from_SP(IOS_ctd.sal);
IOS_ctd.ox=IOS_ctd.ox*44.66;

% Temp prior to Feb 97 in IPTS68
idx=IOS_ctd.mtime<datenum(1997,02,21);
IOS_ctd.temp(:,idx)=gsw_t90_from_t68(IOS_ctd.temp(:,idx));

save('/Users/samuelstevens/Dropbox/Hakai/ctd/data/hecate/ios_ctd.mat',...
    'IOS_ctd');