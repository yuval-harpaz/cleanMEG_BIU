function allCleaned = cleanAllByCoefft(MEG,REF, allCoefs, eigVec, samplingRate, bands)
% find cleaning coeficients and clean all MEG data
% [allCleaned, allCoefs] = coefsAllByFFT(MEG,REF,samplingRate, bands);
%
% MEG          - data of many (all) MEG channel
% REF          - the reference sensors
% samplingRate - in samples per sec.
% bands        - frequency bands to work on [default 0-4,4-8,8-16,...]
%              WARNING the default is based on sampling rate of ~2000/s
% allCleaned   - the cleaned MEG channels (in the time domain)
% allCoefs     - the coeficients for channel and each band in the frequency domain
% eigVec      -  eigen vectors for the REF channels


% Nov-2008 MA

%% initialize
[numMeg,numDataMeg] = size(MEG);
[numRef,numDataRef] = size(REF);
if numDataMeg ~= numDataRef
    error('MEGanalysis:wrongDim', 'data length of REF and MEG must be equal')
end
numData = numDataMeg;

dF = samplingRate/numData;
if ~exist('bands', 'var'), bands=[]; end
if isempty(bands)
%     bands= round([1,4-dF; 4,8-dF; 8,16-dF; 16,32-dF; 32,64-dF; 64,128-dF;...
%         128,256-dF; 256,512-dF; 512,0.9*samplingRate/2-dF]/dF);  
    bands1= round([dF,1-dF; 1,2-dF; 2,4-dF; 4,8-dF; 8,16-dF; 16,32-dF]/dF);
    bb1 = 32:16:0.9*samplingRate/2-16;
    bb2 = (bb1(2):16:0.9*samplingRate/2) -dF;
    bands2 = round([bb1',bb2']/dF);
    bands = [bands1; bands2];
end
% test if enough data in each band
while min(bands(:,2)-bands(:,1))<2*numRef
    b1 = bands(1,:);
    bands(2,1) = b1(1);
    bands(1,:)=[];
    if size(bands,1)<3
        warning('MEGanalysis:minData','notEnough data for cleaning by bands')
    end
end

allCleaned=zeros(size(MEG));
% convet the REFs to the frequency domain
REF = REF-repmat(mean(REF),numRef,1);
% The ref's are not linearly independent project them on principal
% components, use all except the least significant 1
REF =REF'*eigVec;
REF=REF(:,2:end)';
numRef = size(REF,1);
REFfft = zeros(size(REF));
for ii = 1:numRef
    REFfft(ii,:) = fft(REF(ii,:));
end

%% clean all messages
for ii=1:numMeg
    cleanA1 = clean1ByCoefft(MEG(ii,:), REFfft, allCoefs(:,:,ii), eigVec, samplingRate, bands);
    allCleaned(ii,:) = cleanA1;
end

return
