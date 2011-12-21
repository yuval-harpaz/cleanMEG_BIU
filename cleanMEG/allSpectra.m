function [AllPSD, F] = allSpectra(signal,samplingRate,fResolution, method)
% compute power spectra for all channels in signal


%Oct 2008  MA

%% initialize
numChannels = size(signal,1);
[powerSpct, F] = myPSD (signal(1,:), samplingRate, fResolution, method, 'noDC');
numFrequencies = length(F);
AllPSD = zeros(numChannels, numFrequencies);

%% compute
for ii = 1:numChannels
    AllPSD(ii,:) = myPSD(signal(ii,:), samplingRate, fResolution, method, 'noDC');
end
return