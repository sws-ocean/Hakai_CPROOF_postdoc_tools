function ctd=ios_rd(filnam,varargin)
% IOS_RD reads IOS-format CTD files
%  A=IOS_RD(FILENAM) reads a file
%
% 
%  R Pawlowicz 16/Oct/04
%
%  Many updates, latest 13/Apr/2011
%  July/2013 - added E, I format indicators
%  Dec/2018  - some reformatting, other small changes
%  July/2022  - made it work for current meter data

fid=fopen(filnam,'r');

  
zone=nan;

long=NaN;
stn=NaN;
np=NaN;
nchan=NaN;
ofs=27;
formstring='';
padstring=NaN;

% Skip first 2 lines
l=fgetl(fid);l=fgetl(fid);l=fgetl(fid);
 
[fpath,fname,fext]=fileparts(filnam);

ctd=struct('type','IOS CTD','file',[fname fext],'time','','mtime',NaN,...
           'latitude',NaN,'longitude',NaN,'depth',NaN,...
           'station','');
            
 
l=[fgetl(fid) '     '];
while ~strcmp(l(1:4),'*END')
%%     disp(['while->' l]);
  if l(1)=='*'
    switch l(2:5)
        case 'FILE'
           l=[fgetl(fid) '     '];
           while l(1)~='*'
             switch l(5:9)
               case 'START'
                  switch l(ofs+[0:2]),     
                    case {'GMT','UTC','utc'}  
                      zone=0;
                    case 'PST'
                      zone=8;
                    case 'PDT'
                      zone=7;    
                    otherwise
                      zone=0;
                      error(['Unknown time zone ' l(ofs+[0:2])]);
                  end
                  ctd.time=l(ofs+[0:22]);
                          % Not sure what format this reads:  
                  %  ctd.mtime=datenum(l(ofs+[9:13 8 4:7]))+datenum(l(ofs+[15:22]))+zone/24;
                          % Should be   UTC 2007/05/27 21:17:14
                  ctd.mtime=datenum(l(ofs+[4:22]),'yyyy/mm/dd HH:MM:SS')+zone/24;
                 case 'TIME '
                     switch l(10:13)
                         case 'INCR'
                            vals=sscanf(l(ofs:end),'%f');
                            ctd.timeincrement=vals(1)+vals(2)/24+vals(3)/1440+vals(4)/86400 + vals(5)/86400000;
                     end
               case 'NUMBE'
                    switch l(15:21)
                      case 'RECORDS'
                         np=sscanf(l(ofs:end),'%d');
                      case 'CHANNEL'
                         nchan=sscanf(l(ofs:end),'%d');
                    end
               case 'DATA '
                    switch l(10:13)
                        case 'DESC'
                          ctd.type=['IOS ' deblank(l(ofs:end))];
                    end
               case 'FORMA'
                  formstring=upper(l(ofs:end));
                  kk=find(formstring=='(');      % Extract the format
                  ll=find(formstring==')');
                  formstring=formstring(kk+1:ll-1);
               case 'PAD  '
                  padstring=sscanf(l(ofs:end),'%f');      
               case '$TABL'
                  switch l(13:20)
                    case 'CHANNELS'
                          l=fgetl(fid);
                          l=fgetl(fid);
                          % Find starting columns
                          l(1:6)='-';
                          colnums=find(l==' ');
                          span=zeros(2,nchan);
                          for i=1:nchan
                             l=[fgetl(fid) repmat(' ',1,50)];
                             cnam=l((colnums(1)+1):colnums(2));
                             if cnam(1)==''''
                                 cnam=cnam(2:(find(cnam(2:end)=='''')));
                             else
                                 cnam=cnam(1:min(find(cnam(2:end)==' ')));
                             end         
                                    
                          switch lower(cnam)
                                  case {'depth','depth:nominal'}
                                      ctd.varlabel{i}='depth';
                                      pressureindex=i;
                                   case {'pressure','pressure:ctd','pressure:straingauge'}
                                       ctd.varlabel{i}='pr';
                                       pressureindex=i;
                                   case {'temperature','temperature:primary','temperature:reversing',...
                                         'temperature:ctd','temperature:upcast','temperature:high_res'}
                                      ctd.varlabel{i}='temp';  
                                   case {'temperature:secondary'}
                                      ctd.varlabel{i}='temp2';  
                                   case {'temperature:draw'}
                                      ctd.varlabel{i}='temp_draw';  
                                   case {'conductivity_ratio','conductivity','ms/cm','conductivity:primary',...
                                           'conductivity:upcast'}
                                      ctd.varlabel{i}='cond';    
                                   case 'conductivity:secondary'
                                      ctd.varlabel{i}='cond2';   
                                   case {'salinity','salinity:t0:c0','salinity: practical','salinity:bottle',...
                                         'salinity: pre-1978','salinity:t1:c1','salinity:ctd','salinity:practical:t1:c1',...
                                         'salinity:upcast'}
                                      ctd.varlabel{i}='sal'; 
                                   case {'sigma-t','sigma-t:ctd','sigmat-t'}
                                      ctd.varlabel{i}='sigmat';
                                   case {'sigma-theta'}
                                      ctd.varlabel{i}='sigmatheta';
                                   case {'fluorescence','fluorescence:uru','fluorescence:uru:wetlabs',...
                                         'fluorescence:seatech','fluorescence:seapoint','fluorescence:uru:seapoint',...
                                     'fluorescence:ctd','fluorescence:wetlabs','fluorescence:uru:seatech',...
                                     'fluorescence:wetlabs:eco-afl','fluorescence:uru:upcast'}
                                      ctd.varlabel{i}='fluor';
                                   case {'transmissivity','transmissivity:ctd','turbidity:wetlabs'}
                                      ctd.varlabel{i}='xmiss';
                                   case {'transmissivity2'}
                                      ctd.varlabel{i}='xmiss2';
                                   case  {'transmissivity:green'}
                                      ctd.varlabel{i}='xmissGreen';
                                   case {'par','par1'}
                                      ctd.varlabel{i}='par';                        
                                   case {'oxygen','oxygen:dissolved','dissolved_oxygen'}
                                      ctd.varlabel{i}='ox';
                                   case {'oxygen:saturation' }
                                      ctd.varlabel{i}='oxSAT';
                                   case {'oxygen:dissolved:sbe','oxygen:dissolved:ctd'}
                                      ctd.varlabel{i}='oxML';     
                                   case {'phosphate','phosphate(inorg)'}
                                      ctd.varlabel{i}='po4';
                                   case 'nitrite'
                                      ctd.varlabel{i}='no2';
                                   case {'nitrate','nitrate_plus_nitrite'}
                                      ctd.varlabel{i}='no3';
                                   case 'ammonium'
                                      ctd.varlabel{i}='ammonium';         
                                   case {'silicate','silicate:corrected'}
                                      ctd.varlabel{i}='si';
                                   case 'sulphate so4'
                                      ctd.varlabel{i}='so4';      
                                   case {'chlorophyll','chlorophyll:extracted'}
                                      ctd.varlabel{i}='chl';
                                   case {'chlorophyll:extracted:>0.7um'}
                                      ctd.varlabel{i}='chl007';
                                   case {'chlorophyll:extracted:>5.0um'}
                                      ctd.varlabel{i}='chl050';
                                   case 'phaeo-pigment:extracted'
                                      ctd.varlabel{i}='phaeopig';         
                                   case 'production:primary'
                                      ctd.varlabel{i}='production';       
                                   case 'adenosine_triphosphate'
                                      ctd.varlabel{i}='atp';
                                   case {'ph','ph:nbs','ph:sbe','ph:sbe:nominal'}
                                      ctd.varlabel{i}='ph';
                                   case {'alkalinity','alkalinity:total','alkalinity:carbonate'}
                                      ctd.varlabel{i}='ta';      
                                   case 'carbon:dissolved:inorganic'
                                      ctd.varlabel{i}='DIC';                                   
                                   case 'dimethyl_sulphide'
                                      ctd.varlabel{i}='DMS';           
                                   case 'quality_flag:temp'
                                      ctd.varlabel{i}='temp_flag';        
                                   case 'quality_flag:sali'
                                      ctd.varlabel{i}='sal_flag';        
                                   case {'sound velocity','speed:sound','speed:sound:1'}
                                      ctd.varlabel{i}='ssp';
                                   case 'speed'
                                      ctd.varlabel{i}='speed';
                                   case 'direction:geog(to)'
                                      ctd.varlabel{i}='direc_to';
                                   case 'speed:east'
                                      ctd.varlabel{i}='u';
                                   case 'speed:north'
                                      ctd.varlabel{i}='v';
                                   case 'speed:up'
                                      ctd.varlabel{i}='w';
                                  case 'amplitude:beam1'
                                      ctd.varlabel{i}='amp1';
                                  case 'amplitude:beam2'
                                      ctd.varlabel{i}='amp2';
                                  case 'amplitude:beam3'
                                      ctd.varlabel{i}='amp3';
                                   case 'heading'
                                      ctd.varlabel{i}='heading';
                                   case 'pitch'
                                      ctd.varlabel{i}='pitch';
                                   case 'roll'
                                      ctd.varlabel{i}='roll';
                                   case 'quality_flag:soun'
                                      ctd.varlabel{i}='ssp_flag';        
                                   case {'record #','record_number'}
                                      ctd.varlabel{i}='recnum';        
                                   case 'date'
                                      ctd.varlabel{i}='recdate';        
                                   case 'time'
                                      ctd.varlabel{i}='rectime';        
                                    case {'par','par:reference'}
                                      ctd.varlabel{i}='PAR';        
                                    case 'number_of_bin_records'
                                      ctd.varlabel{i}='numbinrecs';        
                                     case 'flag'
                                      ctd.varlabel{i}='flag';        
                                     case 'sample_number'
                                      ctd.varlabel{i}='samplenum';        
                                     case 'bottle:position'
                                      ctd.varlabel{i}='bottlepos';        
                                     case 'bottle_number'
                                      ctd.varlabel{i}='bottlenum';        
                                     case 'bottle:firing_sequence'
                                      ctd.varlabel{i}='bottlesequence';   
                                     case 'reference'
                                      ctd.varlabel{i}='reference';
                                     case 'latitude'
                                      ctd.varlabel{i}='latitude';
                                     case 'longitude'
                                      ctd.varlabel{i}='longitude';
                                     case 'flag:at_sea'
                                      ctd.varlabel{i}='atSea';
                                      atSeaindex=i;
                                     otherwise
                                      kk=[findstr(cnam,':') length(cnam)+1];
                                      ctd.varlabel{i}=lower(cnam(1:kk(1)-1));
                                      if ~strcmp(lower(cnam(1:kk(1)-1)),'flag')
                                          fprintf('Unrecognized channel name: ->%s<-\n',cnam); 
                                      end
                   
                             end
                              ctd.units{i}=deblank(l((colnums(2)+1):colnums(3)));
                              try, span(:,i)=sscanf(l(colnums(3):end),'%f',2);
                              catch span(:,i)=[-Inf Inf]'; end 
                          end
                          l=fgetl(fid);
                    case 'CHANNEL '
                          l=fgetl(fid);l=fgetl(fid);
                          for i=1:nchan
                             l=fgetl(fid);
                             kk=findstr(''' ''',l);
                             if any(kk), for ii=1:length(kk), l(kk(ii)+[0:2])='nan'; end;end
                             kk=findstr('  M',l);
                             if any(kk), for ii=1:length(kk), l(kk(ii)+[0:2])='nan'; end;end
                             kk=findstr('n/a',l);
                             if any(kk), for ii=1:length(kk), l(kk(ii)+[0:2])='nan'; end;end
                             vals=sscanf(l,'%f',4);
                             try, ctd.baddata(i)=vals(2); catch  ctd.baddata(i)=NaN; end
                             if isnan(ctd.baddata(i)), ctd.baddata(i)=padstring; end
                             start(i)=vals(3);
                             width(i)=vals(4);    
                             if isnan(width(i))
                               ll=findstr('F',upper(l));
                               if any(ll)
                                 width(i)=sscanf(l(ll+1:end),'%d'); 
                               end
                               if any(findstr('YYYY/MM/DD',upper(l)))
                                   width(i)=11;
                               elseif any(findstr('HH:MM:SS',upper(l)))
                                   width(i)=9;
                               end
                             end
                          end
                          if isnan(start(1)), start=cumsum([1 width]); end
                          l=fgetl(fid);      
                      otherwise
                          while ~strcmp(l(5:8),'$END')
                             l=fgetl(fid);
                          end
                  end
 
               case '$REMA'
                          while ~strcmp(l(5:8),'$END')
                             l=[fgetl(fid) '         '];
                          end
 
               otherwise
%%                  l=[fgetl(fid) '    '];
             end
             l=[fgetl(fid) '              '];
           end
        case 'LOCA'
           l=[fgetl(fid) '          '];
           while l(1)~='*'
              switch l(5:8)
                case 'LATI'
                      ctd.latitude=sum(sscanf(l(ofs:end),'%f',2)./[1 60]');
                case 'LONG'
                      ctd.longitude=-sum(sscanf(l(ofs:end),'%f',2)./[1 60]');
                case 'WATE'
                      ctd.depth= sscanf(l(ofs:end),'%f',1);
                case 'STAT'
                      ctd.station= deblank(l(ofs:end)); 
                  otherwise
              end
              l=[fgetl(fid) '         '];
           end
 
        case {'ADMI','INST','COMM','HIST','CALI'}
            l=[fgetl(fid) '      '];
            while l(1)~='*'
               l=[fgetl(fid) '     '];
            end
 
        case 'END '
        otherwise
          l=[fgetl(fid) '          '];
    end
     else
      l=[fgetl(fid) '        '];
  end
end

% It is possible that we have multiple columns of the same data
% (i.e. bottle and CTD oxygen), so I'm going to rename the
% labels slightly.

if length(unique(ctd.varlabel))~=length(ctd.varlabel)
  for k=1:length(ctd.varlabel)-1
    kk=strmatch(ctd.varlabel(k),ctd.varlabel(k+1:end));
    if any(kk)
      for l=1:length(kk)
         ctd.varlabel{k+kk(l)}=[ctd.varlabel{k+kk(l)} sprintf('%d',l+1)];
      end
    end
  end
end

% Read the data. Unfortunately the baddata flag sometimes
% removes spaces between numbers so we have to do this the slow way -
% read a line, then take the given format and try to read things in
% that way. Fraught with possibilities of error. All kinds of things have
% been found. Sometimes the # of records does not match with
% that given in header. Sometimes the bad flags overflow the given
% column block. Sometimes the bad data flags are inappropriate 
% (I have noticed that 9999.99 is sometimes written as 9999.9902).
%

data=nan+zeros(np,nchan);
for k=1:np
  l=fgetl(fid);
  if length(l)>1  % Reached end-of-file prematurely
    if isempty(formstring)
     if length(start)==nchan, start=[start length(l)]; end
     for ll=1:nchan
         l=[l repmat(' ',1,10)];
         num=l(start(ll):start(ll+1)-1);
         
         num(num==' ')=[];  
         if any(num=='/')
             data(k,ll)=datenum(num,'yyyy/mm/dd');
         elseif any(num==':')
             data(k,ll)=datenum(num,'HH:MM:SS');
         else
             try
                 data(k,ll)=sscanf(num,'%f'); 
             catch 
                 data(k,ll)=NaN; 
             end    
             if data(k,ll)==9999.9902, data(k,ll)=NaN; end
         end
  
     end
   else
    vals=fortscan(l,formstring);
    ib=vals==9999.9902;
    if any(ib), vals(ib)=NaN; end
    data(k,:)=vals;
    end 
  else
   data(k,:)=nan+zeros(1,nchan);
  end
end
fclose(fid);

% Get rid of lines with no depth data.
% Only can do this if I recognized a pressure channel!
if exist('pressureindex','var')
  data(isnan(data(:,pressureindex)),:)=[];
end

% Now copy the data over to the data structure. Note that
% constant data (i.e. with max==min) will be deleted...
%  July/2022 - made an exception for a pressureindex as sometimes this
%              is set to constant for current meters.
%  Sep/2022 - made another exception for atSea flags.
%
% also data values outside the stated data span will be set to NaN
% under the assumption that it's probably incorrectly decoded BUT
% you might not want this!
%

for k=1:length(ctd.varlabel)
    if isfield(ctd,'baddata')
       data(data(:,k)==ctd.baddata(k),k)=NaN;
    end
    data(data(:,k)<span(1,k) | data(:,k)>span(2,k) ) = NaN;
    if span(1,k)~=span(2,k) | np==1 | (exist('pressureindex','var') && pressureindex==k) ...
            | (exist('atSeaindex','var') && atSeaindex==k)
      ctd.(ctd.varlabel{k})=data(:,k);  
    end
end


function vals=fortscan(l,format)
% Fortran format string

format=upper(format);
comma=[0,find(format==','),length(format)];
n=1;
while n<length(comma)
  switch format(comma(n)+1)
    case 'F'
      fspec=sscanf(format(comma(n)+2:comma(n+1)),'%d.%d');
      vals(n)=sscanf(l(1:fspec(1)),'%f');
      l=l(fspec(1)+1:end);
    case 'E'
      fspec=sscanf(format(comma(n)+2:comma(n+1)),'%d.%d');
      vals(n)=sscanf(l(1:fspec(1)),'%e');
      l=l(fspec(1)+1:end);
    case 'I'
      fspec=sscanf(format(comma(n)+2:comma(n+1)),'%d');
      vals(n)=sscanf(l(1:fspec(1)),'%d');
      l=l(fspec(1)+1:end);
    otherwise
      error(['Unrecognized format specifier: ' format(comma(n)+1:comma(n+1)-1)]);
  end
  n=n+1;
end
       
      

