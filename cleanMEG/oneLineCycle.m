function [meanLine,noiseSamp,cycCount] = oneLineCycle(dataA, whereUp,noiseThr)
% average all cycles of 50Hz
% meanLine = oneLineCycle(dataA, whereUp);
% dataA    - nChannelsXmSamples
% whereUp  - list of indices where the 50Hz goes up
% noiseThr - noise Threshold in z scores

% meanLine - the mean triggered on the 50hz up stroke
% noiseSamp- are samples which are members of noisy cycles
% cycCount - the number of good cycles


% Sep-2008  MA
% May-2014  YH added noise rejection
%% initialize
if ~exist('noiseThr','var')
    noiseThr=[];
end
if isempty(noiseThr)
    noiseThr=5;
end
    
maxL = max(diff(whereUp));
[nChannels,mSamples] = size(dataA);
numCycles = floor(mSamples/maxL)-1;
% numCycles=1000;
sum1 = zeros(nChannels,maxL+1);
% Estimate noise
for cycle = 1:(numCycles-2)
    startCycle = whereUp(cycle);
    amp1(cycle) = mean(abs(dataA(startCycle:startCycle+maxL)-mean(dataA(startCycle:startCycle+maxL)))); %#ok<AGROW>
end
amp2=(amp1-mean(amp1))./std(amp1);
noise=min(amp1(amp2>=noiseThr));
if isempty(noise)
    noise=max(amp1); % to accept all segments
end
% average
cycCount=0;
noiseSamp=[];
for cycle = 1:numCycles
    startCycle = whereUp(cycle);
    if startCycle+maxL <= size(dataA,2)
        if mean(abs(dataA(startCycle:startCycle+maxL)-mean(dataA(startCycle:startCycle+maxL))))<=noise
            cycCount=cycCount+1;
            sum1 = sum1 + dataA(:,startCycle:startCycle+maxL);
        else
            noiseSamp=[noiseSamp,startCycle:(startCycle+maxL)]; %#ok<*AGROW>
        end
    end
end

meanLine = sum1/cycCount;
