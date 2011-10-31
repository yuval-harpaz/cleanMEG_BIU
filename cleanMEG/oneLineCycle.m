function meanLine = oneLineCycle(dataA, whereUp)
% average all cycles of 50Hz
% meanLine = oneLineCycle(dataA, whereUp);
% dataA    - nChannelsXmSamples
% whereUp  - list of indices where the 50Hz goes up
% meanLine - the mean triggered on the 50hz up stroke

% Sep-2008  MA

%% initialize
maxL = max(diff(whereUp));
[nChannels,mSamples] = size(dataA);
numCycles = floor(mSamples/maxL)-1;
% numCycles=1000;
sum = zeros(nChannels,maxL+1);

for cycle = 1:numCycles
    startCycle = whereUp(cycle);
    if startCycle+maxL <= size(dataA,2)
        sum = sum + dataA(:,startCycle:startCycle+maxL);
    end
end

meanLine = sum/numCycles;

