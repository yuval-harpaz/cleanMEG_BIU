function  [xClean, mean1] = cleanWphase(x,whereUp,interpNo)
% clean periodic signal from x, allow phase precession


% Oct-2010  MA

%% initialize
if ~exist('interpNo','var'), interpNo=[]; end
if isempty(interpNo), interpNo = 10; end

xInt = interp(x,interpNo);
whereUpInt = interpNo*(whereUp-1) +1;


%% find the mean
[phasePerSample, phaseAtTrig] = findPhasePrecession(whereUpInt);
[mean3cycles, zeroAt] = findMeanCycle(xInt, whereUpInt, 1);

%% clean 1 cycle at a time
for iTrig = 1:length(whereUp)-1
    xI0 = whereUpInt(iTrig);
    xI1 = whereUpInt(iTrig+1)-1;
    phase0 = phaseAtTrig(iTrig);
    % phase1 = phaseAtTrig(iTrig+1);
    % find the closest index
    i0 = round(phase0/phasePerSample)+zeroAt;
    i1 = i0 +(xI1-xI0);
    %i1 = round((phase1-zeroAt)/phasePerSample);
    xInt(xI0:xI1) = xInt(xI0:xI1)-mean3cycles(i0:i1);
end


%% treat the edges
% the beginning
if whereUpInt(1)>1
    xI0 = 1;
    xI1 = whereUpInt(1)-1;
    phase0 = phaseAtTrig(1)- xI1*phasePerSample;
    % phase1 = phaseAtTrig(1);
    i0 = round(phase0/phasePerSample) +zeroAt;
    i1 = i0 +(xI1-xI0);
    xInt(xI0:xI1) = xInt(xI0:xI1)-mean3cycles(i0:i1);
end
% the end
if whereUpInt(end)<length(xInt)
    xI1 = length(xInt);
    xI0 = whereUpInt(end)+1;
    % phase1 = xI1*phasePerSample -phaseAtTrig(end);
    phase0 = phaseAtTrig(end);
    i0 = round(phase0/phasePerSample) +zeroAt;
    i1 = i0 +(xI1-xI0);
    xInt(xI0:xI1) = xInt(xI0:xI1)-mean3cycles(i0:i1);
end

%% wrap up
xClean = xInt(1:10:end);
mean1 = mean3cycles;

return
