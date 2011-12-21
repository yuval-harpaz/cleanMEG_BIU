function [noHiFnoise, excess] = findOKperiods(x, samplingRate, noiseLimit, tStep, lineFreq, toPlot)
%  Find 10 sec periods without big spectral lines at hi frequency
%    noHiFnoise = findOKperiods(x, samplingRate, noiseLimit, tStep,...
%    lineFreq);
%
% x            - raw vector with MEG data
% samplingRate - in samples per second
% noiseLimit   - max allowed ratio of hi-f peaks relative to mean of lower
%                frequencies [default 180]
% tStep        - steps to advance the analysis, in seconds [default 5 sec]
% lineFreq     - [default 50]
% toPlot       - 'PLOT' if to plot a histogram of the outliers values
%                [default 'NoPlot']
%
% noHiFnoise - 2xN matrix with start and end points of periods with no high
%              frequency noise.
% excess     - the ratio of peak in hi frequencies to MAD of lower
%              frequencies
%
%   ALGORITHM
% The PSD of pieces 10 s long is computed.  The median and mad of
% frequencies between 103 and 147 Hz is computed.  If the highest peak for
% frequencies above 160 Hz is above median+noiseLimit*mad the piece is
% considered noisy.  the 10 s pieces are advanced in tStep and the
% computation repeats.

% Dec 2011  MA

%% internal params
tPiece = 10;         % analyze pieces in length of 10 sec
% tStep = 5;           % move time at steps of 5 sec
fResolution = 0.25;  % resolution of power spectra
% noiseLimit = 100;    % how many MADs is the peak at hi freq

%% initialize
if samplingRate<600
    error('MATLAB:MEGanalysis:badParam',...
        ['Sampling rate must be greater then 700 to avoid ' ...
        'confusion with line frequency artifacts'])
end
if ~exist('tStep', 'var'), tStep=[]; end
if isempty(tStep), tStep=tPiece/2; end
if ~exist('lineFreq', 'var'), lineFreq=50; end
if isempty(lineFreq), lineFreq=50; end
if ~exist('noiseLimit', 'var'), noiseLimit=[]; end
if isempty(noiseLimit), noiseLimit=180; end
if ~exist('toPlot', 'var'), toPlot=[]; end
if isempty(toPlot), toPlot='noPlot'; end
toPlot = strcmpi(toPlot, 'PLOT');
lf = lineFreq;
maxF = ceil(0.7*samplingRate/2);
% frequency range for base line
F1=2*lf+3;
F2=3*lf-3;
% hi frequencies excluding lf harmonics
FF = 3*lf:lf:0.9*maxF;
Fstrt = FF+10;
Fends = [FF(2:end)-10 maxF];
if (Fends(end)-Fstrt(end)) <30  % ignore the last piece
    Fstrt(end)=[];
    Fends(end)=[];
end
numHiHarmonics = length(Fstrt);

numData = size(x,2);

tTotal= numData/samplingRate;
numSteps = floor((tTotal-tPiece)/tStep);
hiFnoise = false(1, numSteps);
excess = nan(1, numSteps);

%% find frequency indices to analyze
t0 = 0;
t1 = tPiece;
i0 = floor(t0*samplingRate)+1;
i1 = ceil(t1*samplingRate);
y = x(i0:i1);
[~, F] = myPSD (y, samplingRate, fResolution, 'FFT', 'noDC');
% baselineFreq
f1 = find(F>F1,1);
f2 = find(F>F2,1);
% hi frequencies
f3 = nan(1,numHiHarmonics);
f4 = f3;
for ii = 1:numHiHarmonics
    f3(ii) = find(F>Fstrt(ii),1);
    f4(ii) = find(F>Fends(ii),1);
end
hiFindx = false(1,length(F));
for ii = 1:numHiHarmonics
    hiFindx(f3(ii):f4(ii)) = true;
end

%% search for hiFreq artifacts
for ii = 1:numSteps
    t0 = (ii-1)*tStep;
    t1 = t0+ tPiece;
    if ii == numSteps  % the last piece include all
        t1 = tTotal;
    end
    i0 = floor(t0*samplingRate)+1;
    i1 = ceil(t1*samplingRate);
    y = x(i0:i1);
    
    %% try to see if there is high frequncy noise in this piece
    PSD = myPSD (y, samplingRate, fResolution, 'FFT', 'noDC');
    mx=max(PSD(hiFindx));  %max except lf harmonics
    mn=median(PSD(f1:f2));
    ss=mad(PSD(f1:f2));
    oPSD =(mx-mn)/ss;
    excess(ii) = oPSD;
end

%% mark pieces with high noise
for ii = 1:numSteps
    oPSD = excess(ii);
    if oPSD < noiseLimit
        hiFnoise(ii) = false;
    else
        hiFnoise(ii) = true;
    end
end

%% summerize and find if any noisy
thisChan = hiFnoise;
if any(thisChan)  % any noisy peice here
    [runLength, nameId] = findRuns(thisChan);
    lrl = length(runLength);
    if lrl==1  % all data in the same type
        if nameId==0  % all OK
            noHiFnoise = [0,tTotal];
        else          % all noisy
            noHiFnoise = [];
        end
        return
    end
    strtIndx = [0 cumsum(runLength)];
    quietStrt = nan(1,lrl);
    quietEnd = quietStrt;
    % does it start as a quiet piece?

    if nameId(1)==0
        quietStrt(1) = 0;
        quietEnd(1) = tStep*(strtIndx(2)-1) +tPiece;
        kk =3;
        jj =2;
        while kk<length(strtIndx)-1
            quietStrt(jj) = tStep*(strtIndx(kk));
            quietEnd(jj) = tStep*(strtIndx(kk+1)-1) +tPiece;
            kk = kk+2;
            jj = jj+1;
        end
        % treat the last one
        if (nameId(end)==0)&&(kk<=lrl)  % last one is also quiet
            quietStrt(jj) = tStep*strtIndx(end-1);
            quietEnd(jj) = tTotal;
        end
    else  % starts with noise
        quietStrt(1) = tStep*(strtIndx(2));
        quietEnd(1) = tStep*(strtIndx(3)-1) +tPiece;
        kk = 4;
        jj = 2;
        while kk<length(strtIndx)-1
            quietStrt(jj) = tStep*(strtIndx(kk));
            quietEnd(jj) = tStep*(strtIndx(kk+1)-1) +tPiece;
            kk = kk+2;
            jj = jj+1;
        end
        % treat the last one
        if (nameId(end)==0)&&(kk<=lrl)  % last one is also quiet
            quietStrt(jj) = tStep*strtIndx(end-1);
            quietEnd(jj) = tTotal;
        end
    end
    quietStrt(isnan(quietStrt)) = [];
    quietEnd(isnan(quietEnd)) = [];
    % merge overlapping regions
    I = find(quietStrt(2:end)<=quietEnd(1:end-1));
    if ~isempty(I)
        quietEnd(I) = [];
        quietStrt(I+1) = [];
    end
    noHiFnoise = [quietStrt ; quietEnd];
else % all was quiet
    noHiFnoise = [0; tTotal];
end

%% plot histogram if needed
if toPlot
    figure
    if max(excess) < noiseLimit
        hist(excess,30);
    else
        hist(excess, 4*max(excess)/noiseLimit);
    end
end

return
