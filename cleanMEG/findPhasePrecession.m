function [phasePerSample, phaseAtTrig] = findPhasePrecession(whereUp)
% find precession between LF and samples.
% 2-nd trig is considered as 0 phase
%  [phasePerSample, phaseAtTrig] = findPhasePrecession(whereUp);
%
% whereUp        - list of sample numbers where the trigger bit flipped up
%
% phasePerSample - how many radians per one sample
% phaseAtTrig    - what is the phase at each trig, in the range of +/-pi

% Oct-2010  MA

%% initialize
numCycles = length(whereUp) -1;
samplesPerCycle = (whereUp(end)-whereUp(1))/numCycles;
phasePerSample = 2*pi/samplesPerCycle;

%% find the phase per trig in the range of -pi to pi
% the phase of the 3-rd trig is considered as 0
phaseAtTrig = (whereUp-whereUp(2))*phasePerSample;
% bring to +/-pi
phaseAtTrig = mod(phaseAtTrig,2*pi);
big = phaseAtTrig>=pi;
phaseAtTrig(big) = phaseAtTrig(big)-2*pi;

return