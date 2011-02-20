function [where, amplitude, atEnd] = findBigStep(Data, samplingRate, outLier, stepDur, maxAmplitude)
%  Find where is a big step in Data
%   [where, amplitude] = findBigStep(Data,samplingRate,outLier,stepDur, maxAmplitude);
%
% Data         - vector with data
% samplingRate - in samples per sec.
% outLier      - how many MADs is considered an outlier [default 10]
% stepDur      - what fraction of samplingRate data should be above peak/3
%                to be considered as step.  [default 1/5]
% maxAmplitude - amplitude at step peak at least so big. [default 1e-11]
% where        - index of the peak immediately after the step (empty if no
%                step found).
% amplitude    - amplitude of the peak immediately after the step (empty 
%                if no step found).
% atEnd        - the step is too close to the end

% Dec-2010  MA

%% initialize
if ~exist('outLier', 'var'), outLier =[]; end
if isempty(outLier), outLier = 10; end
if ~exist('stepDur', 'var'),stepDur  =[]; end
if isempty(stepDur), stepDur = 1/5; end
if ~exist('maxAmplitude', 'var'), maxAmplitude =[]; end
if isempty(maxAmplitude), maxAmplitude = 1e-11; end
backPieceSize = ceil(2*stepDur*samplingRate);
forePieceSize = ceil(1.5*stepDur*samplingRate);
numData = length(Data);
where=[];
amplitude=[];
atEnd = false;
%  medData = median(Data);

%% try repetitively
while true
    absData = abs(Data);
    dData = diff(Data);
    absDdata = abs(dData);
    medD = median(dData);
    madD = mad(dData);
    if madD==0
        return
    end
    
    %% find biggest jump
    maxDd = max(absDdata);
    if (maxDd-medD)/madD <outLier  % no big jump
%         where = wherePeak;
%         amplitude = Data(where);
        return
    end
    
    whereJump = find(absDdata ==maxDd,1);
    % whereJump = find(absDdata ==maxDd);
    
    %% decide if a step function
    I0 = whereJump-10;
    if I0<1, I0=1; end
    I1 = whereJump+10;
    if I1>numData, I1=numData; end
    K0=I0;
    K1 = K0 + forePieceSize; %~ 5sec for 0.1 Hz
    if K1>numData, K1=numData; end
    localData = Data(K0:K1);
    maxAmp = max(absData(I0:I1));  
    wherePeak = find(abs(localData)==maxAmp,1)+K0 -1;
    % wherePeak = find(abs(localData)==maxAmp)+K0 -1;
    peakAmplitude = absData(wherePeak);
    J0 = wherePeak-backPieceSize;
    if J0<1, J0=1; end
    J1 =  wherePeak-3;
    baseLine = median(Data(J0:J1));
    if all( (abs(Data(wherePeak :K1)-baseLine)) >peakAmplitude/3) ...
            && peakAmplitude >= maxAmplitude  % this is a step
        where = wherePeak;
        amplitude = Data(where);
        if numData-wherePeak<forePieceSize, atEnd=true; end
        return
    else % not a real step search again
        I0 =wherePeak-10;
        if I0<1, I0=1; end
        I1 =wherePeak+forePieceSize;
        if I1>numData, I1=numData; end
        Data(I0:I1)= baseLine;
    end
end

