function [powerSpct, F] = myPSD (x, samplingRate, fResolution, method, noDC)
% compute power density, either by average of many FFTs or by PSD
% x           - the data in a row vector
% samplingRate - In samples per second
% fResolution - requested frequency resolution in Hz (default 1);
% method      - 'PSD' or 'FFT' (default FFT)
% noDC        _ 'noDC' or 'All' (default noDC)
% powerSpct   - the power spectra densities
% F           - the frequencies
%
% NOTES: 
%   1.  When noDC is specified, the mean is subtracted from each piece (in
%       FFT), or from the whole (in PSD).  In addition the zero's and 1st
%       frequencies are set to the mean of frequencies 3:13.
%   2.  In FFT only one side is returned
%   3.  In FFT the power is forced to be as in the data

% Sep 2008  MA
%  UPDATES
%    Nov-2010 logical tests on frequency resolution versus data length  MA

%% Initialize
if ~exist('fResolution', 'var'), fResolution=[]; end
if isempty(fResolution), fResolution= 1; end
if ~exist('method', 'var'), method=[]; end
if isempty(method), method= 'FFT'; end
if ~strcmpi(method,'FFT') && ~strcmpi(method,'PSD')
    error('MEGanalysis:FrequenctMethod','method must be either FFT or PSD');
end
if ~exist('noDC', 'var'), noDC=[]; end
if isempty (noDC), noDC='noDC'; end
if strcmpi(noDC, 'noDC')
    subtractDC=true;
elseif strcmpi(noDC,'All')
    subtractDC=false;
else
    error('MEGanalysis:FrequenctMethod','subtract DC must be either noDC or All');
end

totalSamples = length(x);

%% test for fit of params
totalT = totalSamples/samplingRate;
if 1/totalT > 1.5*fResolution %error cannot
    warning('MATLAB:MEGanalysis:ConflictingParams'...
        ,'Data too short for resolution of %f cycles', fResolution);
    fResolution = 1.5/totalT;
    disp(['Changed to ' num2str(fResolution) ' Hz']);
end
samples2analyze= round(samplingRate*1/fResolution);

%% break to pieces if FFT
if strcmpi(method, 'FFT')
    scalFactor = 1/mean(abs(diff(x)));
    xScaled = x*scalFactor;
    numSamples = totalSamples;
    numSections = ceil(numSamples/samples2analyze);
    numInF = ceil(samples2analyze/2);
    sectionStart = zeros(1,numSections);
    sectionEnds  = sectionStart;
    for ii = 1:numSections
        sectionStart(ii) = round((ii-1)*numSamples/numSections)+1;
        sectionEnds(ii)  = sectionStart(ii) + samples2analyze -1;
    end
    if sectionEnds(end) > totalSamples
        sectionEnds(end) = totalSamples;
        sectionStart(end) = sectionEnds(end) -samples2analyze +1;
    end
    powerSpct = zeros(1 , numInF);
    F = 0:fResolution:samplingRate/2;
    % ?? may be off by 1 ??
    if length(F)<size(powerSpct,2)
        F = [F, F(end)+fResolution];
    elseif length(F)>size(powerSpct,2)
        F = F(1:end-1);
    end
end

%% do the FFT
if strcmpi(method, 'FFT');
    for sectionNo=1:numSections
        I=(sectionStart(sectionNo):sectionEnds(sectionNo));
        y = xScaled(I);
        if subtractDC
            y = y-mean(y);
        end
        pwr = abs(fft(y));
        pwr = pwr.*pwr;
        powerSpct = powerSpct +pwr(1:numInF);
    end
    powerSpct = var(x)*powerSpct/sum(powerSpct(2:end));
    if subtractDC
        powerSpct(1:2) = mean(powerSpct(3:12));
    end

elseif strcmpi(method, 'PSD')
    Hs = spectrum.welch('Hamming',samples2analyze,50);
    Hpsd = psd(Hs,x,'Fs',samplingRate);
    powerSpct = Hpsd.data';
    if subtractDC
        powerSpct(1) = mean(powerSpct(3:7));
    end
    F=Hpsd.Frequencies;
end

