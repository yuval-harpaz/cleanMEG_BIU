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
meanCycleLength = length(meanCycle);

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
if whereCycle(1)>(2*meanCycleLength-zTime)
    error('MATLAB:MEGanalysis:wrongParameters','first cycle too long')
end
if (length(Data)-whereCycle(end))>(2*meanCycleLength-zTime)
    error('MATLAB:MEGanalysis:wrongParameters','last cycle too long')
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


stdIafter = length(meanCycle)-zTime-1;
cData = Data;
longMean = [meanCycle, meanCycle, meanCycle]; % 3 concatenated cycles
longZero = length(meanCycle)+zTime;           % where is zTime in longmean

%% clean the very first cycle
zData = whereCycle(1);
dd1=1;
dd2 = zData+stdIafter;
d1 = longZero-zData+1;
if d1<1 %error
    d1 =1;
    d2 =dd2-dd1+1;
else
    d2 = longZero+stdIafter;
end
cData(dd1:dd2) = Data(dd1:dd2) - longMean(d1:d2)*Amplitude(1);


%% clean all cycles in the middle
for cycle = 2:length(dCycle)
    dur2 = dCycle(cycle);
    if dur2>meanCycleLength % avoid going into next one
        if dur2>2*meanCycleLength-zTime %  a cycle was missed
            unCleaned(unCleanIndx) = whereCycle(cycle-1);
            unCleanIndx =  unCleanIndx+1;
        end
        dur2=meanCycleLength;
    end
    dd1 = double(whereCycle(cycle)-zTime);
    dd2 = double(dd1+dur2);
    d1 = meanCycleLength;
    d2 = double(d1+dur2);
    cData(dd1:dd2) = Data(dd1:dd2) -longMean(d1:d2)*Amplitude(cycle);
end
%% clean the last cycle
where = whereCycle(end);  %BUG fix in MATLAB??
dur2 = length(Data) -where;
if dur2<meanCycleLength  % a very short piece
    dd1 = double(where-zTime);
    dd2 = double(dd1+dur2);
    d1 = meanCycleLength;
    d2 = double(d1+dur2);
else % long time up to end of data
    dur2 = meanCycleLength;
    if dur2>2*meanCycleLength-zTime %  a cycle was missed
        unCleaned(unCleanIndx) = whereCycle(cycle-1);
        unCleanIndx =  unCleanIndx+1;
    end
    dd1 = double(where-zTime);
    dd2 = double(dd1+dur2);
    d1 = meanCycleLength;
    d2 = d1+dur2;
end
cData(dd1:dd2) = Data(dd1:dd2) -longMean(d1:d2)*Amplitude(cycle);
unCleaned(unCleanIndx:end)=[];

return
