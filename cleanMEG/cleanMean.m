function [cData, unCleaned] = cleanMean(Data, whereCycle, meanCycle, zTime, Amplitude)
% clean the meanCycle from Data
%   cData = cleanMean(Data, whereCycle, meanCycle, Ztime, Amplitude);
%
% Data       - a vector with data to be cleaned
% whereCycle - list of indices into Data where the cycles are
% meanCycle  - the shape of the mean to be subtracted from Data
% Ztime      - where in meanCycle is time 0.  [default 1]
% Amplitude  - list of scaling factors for the meanCycle one per cycle.
%              [default (1,1,...)].
%
% cData      - the cleaned data
% unCleaned  - list of indices where could not clean
% Notes
%  whereCycle may not be exactly periodic, however, mean Cycle must be
%  at least as long as the longest interval in whereCycle.
%  Data, meanCycle must both be raw vectors

% Nov-2009  MA
% UPDATES
%  Dec-2009  fitting algorithm changed  MA
%  Jan-2010  Cycles which are shorter then zTime eliminated.  MA

%% fudge for testing
% zTime = zTime-16;

%% initialize
if ~exist('zTime','var'), zTime=[]; end
if isempty(zTime),  zTime=1; end
if ~exist('Amplitude','var'), Amplitude=[]; end
if isempty(Amplitude), Amplitude =ones(size(whereCycle)); end
unCleaned=zeros(1,length(whereCycle));
unCleanIndx=1;

% logical tests
if size(Data,1)~=1
    error('MATLAB:MEGanalysis:wrongParameters','Data must be a raw vector')
end
if size(meanCycle,1)~=1
    error('MATLAB:MEGanalysis:wrongParameters','meanCycle must be a raw vector')
end
if whereCycle(end)>length(Data)
    error('MATLAB:MEGanalysis:wrongParameters','Indices in whereCycle too big')
end
if whereCycle(1)>(2*length(meanCycle)-zTime)
    error('MATLAB:MEGanalysis:wrongParameters','first cycle too big')
end
if (length(Data)-whereCycle(end))>(2*length(meanCycle)-zTime)
    error('MATLAB:MEGanalysis:wrongParameters','last cycle too big')
end
if length(Amplitude)~=length(whereCycle)
    error('MATLAB:MEGanalysis:wrongParameters','Amplitude and whereCycle do not match')
end
dCycle = diff(whereCycle);
% erase cycles which are too short
while min(dCycle)<zTime
    tooSmall = find(dCycle<zTime,1);
    whereCycle(tooSmall)=[];
    Amplitude(tooSmall)=[];
    dCycle = diff(whereCycle);
end


stdIafter = length(meanCycle)-zTime;
cData = Data;
longMean = [meanCycle, meanCycle, meanCycle]; % 3 concatenated cycles
longZero = length(meanCycle)+zTime;           % where is zTime in longmean

%% clean the very first cycle
zData = whereCycle(1);
dd1=1;
dd2 = zData+stdIafter;
d1 = longZero-zData+1;
if d1<1 %error
    
end
d2 = longZero+stdIafter;
cData(dd1:dd2) = Data(dd1:dd2) - longMean(d1:d2)*Amplitude(1);


%% clean all cycles in the middle
for cycle = 2:length(dCycle)
    dur2 = dCycle(cycle);
    zData = whereCycle(cycle);
    dd1=dd2 +1;
    dBack = zData-dd1;
    dFore = dur2-dBack;
    dd2= dd1 +dur2 -1;
    d1 = longZero -dBack;
    d2 = longZero +dFore -1;
    if dd2>length(cData)
        ddd = dd2-length(cData);
        dd2 = dd2-ddd;
        d2  = d1  +dd2-dd1;
    end
    cData(dd1:dd2) = Data(dd1:dd2) -longMean(d1:d2)*Amplitude(cycle);
end
%% clean the last cycle
dEnd = length(Data)-dd2;
if dEnd>0
    dd1 = dd2+1;
    zData = whereCycle(end);
    d1 = longZero - (zData-dd1);
    d2 = d1+dEnd-1;
    cData(dd1:end) = Data(dd1:end) -longMean(d1:d2)*Amplitude(end);
end
% if dEnd<=length(meanCycle) % use only a piece of mean
%     % dur2 = dCycle(end);
%     zData = whereCycle(end);
%     dd1=dd2+1;
%     dBack = zData-dd1;
%     dFore = length(Data)-zData;
%     d1 = zTime -dBack;
%     d2 = zTime +dFore;
%     meanHere=Amplitude(end)*meanCycle(d1:d2);
%     cData(dd1:end) = Data(dd1:end)-meanHere;
% else  % stitch the begining to the end
%     tmpMean = [meanCycle, meanCycle(1:dEnd-length(meanCycle))];
%     zData = whereCycle(end);
%     dd1=dd2+1;
%     dBack = zData-dd1;
%     d1 = zTime-dBack;
%     dFore = length(Data)-zData;
%     d2= dBack+dFore+1;
%     if d2>length(tmpMean)
%         % warning('MATLAB:MEganalysis:indexTooBig','at %d Cycle is too big - not cleaned', zData)
%         if unCleanIndx>1
%             if unCleaned(unCleanIndx-1)~=zData
%                 unCleaned(unCleanIndx) = zData;
%                 unCleanIndx= unCleanIndx+1;
%             end
%         else
%             unCleaned(unCleanIndx) = zData;
%             unCleanIndx= unCleanIndx+1;
%         end
%         cData(dd1:end) = Data(dd1:end);
%     else
%         cData(dd1:end) = Data(dd1:end)-Amplitude(end)*tmpMean(d1:d2);
%     end  
% end
unCleaned(unCleanIndx:end)=[];

return
