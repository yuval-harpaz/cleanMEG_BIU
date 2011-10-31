function [mean3cycles, zeroAt] = findMeanCycle(x, whereUp, interpNo)
% compute the mean of 3 cycles allowing for phase precession
%    [mean3cycles, zeroAt] = findMeanCycle(x, whereUp);
%   
%
% x           - array with data
% whereUp     - at which sample approximately did the cycle start
% interpNo    - how much to interpolate [default 10]
% 
% mean3cycles - 3 consecutives mean cycles
% zeroAt      - at which sample is zero phase
%
% NOTES: the data (and means) are interpolated 10 fold.
%        the 3-rd trig in whereUp is considered to be at 0 phase

% Oct-2010  MA

%% initialize
if ~exist('interpNo','var'), interpNo=[]; end
if isempty(interpNo), interpNo = 10; end
if interpNo>1
    x = interp(x,interpNo);
    whereUp = interpNo*(whereUp-1) +1;
end
[phasePerSample, phaseAtTrig] = findPhasePrecession(whereUp);
pi2 = 2*pi;
numBack = floor(pi2/phasePerSample)-1;
numFore = floor(2*pi2/phasePerSample)-1;
numSamples = numBack+numFore +1;
% numSamples = ceil(3*pi2/phasePerSample);
sumV = zeros(1,numSamples);
numAveraged = sumV;

%% compute the soum over 3 cycles
for trigNo = 2:length(whereUp)-3
    trig0 = trigNo-1;
    trig3 = trigNo+2;
    phase0 = phaseAtTrig(trig0)-pi2;
    phase3 = phaseAtTrig(trig3)+2*pi2;
    allPhases = phase0:phasePerSample:phase3;
    zeroSample = find(allPhases>=0,1);
%     % debug section
%     if ~isempty(find(diff(allPhases)<=0,1))
%         disp(['Order reversal at ' num2str(trigNo)])
%     end
    indx = ((zeroSample-numBack):(zeroSample+numFore)) +whereUp(trig0);
    sumV = sumV + x(indx);
    numAveraged = numAveraged+1;
end

%% wrap up
mean3cycles = sumV./numAveraged;
zeroAt = numBack;
return
