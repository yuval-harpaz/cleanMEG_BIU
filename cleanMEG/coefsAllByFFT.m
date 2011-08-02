function [allCleaned, coefStruct] = coefsAllByFFT(MEG,REF,samplingRate, bands, prevCoefs, overFlowTh, rot90)
% find cleaning coeficients and clean all MEG data
% [allCleaned, coefStruct] = coefsAllByFFT(MEG,REF,samplingRate, bands, prevCoefs);
%
% MEG          - data of many (all) MEG channel
% REF          - the reference sensors
% samplingRate - in samples per sec.
% bands        - frequency bands to work on [default 1-4,4-8,8-16,...]
%              WARNING the default is based on sampling rate of ~2000/s
% prevCoefs    - use predefined coefs to subtract the REF from MEG
%                it is a structure with the following fields:
%                prevCoefs -  the matrix of cleaninig coeficients
%                eigVec    - the eigen vectors for the REF channels
%                bands     - frequency bands to work on 
%                            [default 1-4,4-8,8-16,...] (may be empty)
% overFlowTh   - threshold for overflow [default 1e-11];
% rot90        - true if to rotate the REF FFT by 90 degrees to accomodate
%                also phase shifts of acceleration. [default false]
%
% allCleaned   - the cleaned MEG channels (in the time domain)
% coefStruct   - structure with allCoefs, eigVec, bands, as follows:
%     allCoefs     - the coeficients for channel and each band in the frequency domain
%     eigVec       - eigen vectors of ref in columns (Only if ref is in time
%                    domain)  It is scaled up by rFactor.
%     bands        - the bands Used in Hz, inclusive bounderies

% Nov-2008 MA
%  UPDATES
% Dec-2008  Lower band changed from 0-4 Hz into 1-4 Hz.
% Jun-2009  Default bands were wrong for low sampling rates - fixed  MA
% Jul-2009  When very large deflections a cahnnel is not cleaned  MA
% Jul-2010  If no. of REFs smaller then 6 use all of them as is
%           toRot for generating 90 deg phase shift is added  MA

%% initialize
[numMeg,numDataMeg] = size(MEG);
[numRef,numDataRef] = size(REF);
if numDataMeg ~= numDataRef
    % under some rare conditions one of them may have an extra sample
    if numDataMeg > numDataRef
        MEG(:,end)=[];
    elseif numDataMeg < numDataRef
        REF(:,end)=[];
    end
    [numMeg,numDataMeg] = size(MEG);
    [numRef,numDataRef] = size(REF);
    if numDataMeg ~= numDataRef
        error('MEGanalysis:wrongDim', 'data length of REF and MEG must be equal')
    end
end
numData = numDataMeg;

if ~exist('overFlowTh', 'var'), overFlowTh=[]; end
if isempty(overFlowTh), overFlowTh = 1e-11; end

if ~exist('prevCoefs', 'var'), prevCoefs=[]; end
if isempty(prevCoefs)
    findNewCoefs=true;
    eigVec=[];  % so it will be recomputed
else
    allCoefs=prevCoefs.allCoefs;
    eigVec  = prevCoefs.eigVec;
    bands   = prevCoefs.bands;
    % check for matching dimensions:
    if size(allCoefs,1)~=size(bands,1) || ...
            size(allCoefs,2)~=size(eigVec,1)-1 || ...
            size(allCoefs,3)~=numMeg || ... 
            size(eigVec,1) ~= numRef
        warning('MEGanalysis:wrongDim',...
            ['Mismatch between dims of allCoefs, eigVec, bands\n'...
            'Coeficients will be recomputed for this piece.'])
        allCoefs=[];
        eigVec=[];
    end
    if isempty(allCoefs)||isempty(eigVec)
        findNewCoefs=true;
    else
        findNewCoefs=false;
    end
end
if ~exist('rot90', 'var'), rot90=[]; end
if isempty(rot90), rot90 = false; end

dF = samplingRate/numData;
if ~exist('bands', 'var'), bands=[]; end
if isempty(bands)  %|| size(bands,1)<3
    warning('MEGanalysis:missingInput:settingBands','Setting bands to their default values');
    warning('OFF', 'MEGanalysis:missingInput:settingBands')
%     bands= round([1,4-dF; 4,8-dF; 8,16-dF; 16,32-dF; 32,64-dF; 64,128-dF;...
%         128,256-dF; 256,512-dF; 512,0.9*samplingRate/2-dF]/dF);  
    bands1= round([dF,1-dF; 1,2-dF; 2,4-dF; 4,8-dF; 8,16-dF; 16,32-dF]/dF);
    bb1 = 32:16:0.9*samplingRate/2-16;
    bb2 = (bb1(2):16:0.9*samplingRate/2) -dF;
    bands2 = round([bb1',bb2']/dF);
    bands = [bands1; bands2];
    % truncate hi ranges if overlaps
    while bands(end,1)>=bands(end,2)-10
        bands(end,:)=[];
        if ~findNewCoefs  % band structure had ti be changed cannot use old
            % cleaning factors
            findNewCoefs=true;
            warning('MEGanalysis:notEnoughData', 'Cannot use old cleaning factors as bands had to be changed')
        end
    end
    bands(end,2) = round((0.9*samplingRate/2-dF)/dF);
% elseif bands(2,1)==bands(1,2)   % dF has to be subtracted
%     if bands(3,1)<1/dF  % needs to be converted to dF units
%         bands = round(bands/dF);
%     end
%     bands(1:end-1,2) = bands(1:end-1,2)-1;  % we assume here that bands is not in cycles but in dF units
elseif size(bands,2)==2  % both limits are given
    bands = round(bands/dF);
    % check for no overlaps between bands
    for ii = 2:size(bands,1)
        if bands(ii,1)>=bands(ii-1,2)  % there is overlap
            bands(ii-1,2) = bands(ii,1)-1;
        end
    end
elseif size(bands,2)==1  % only limits are given
    tempBands = zeros(size(bands,1)-1,2);
    for ii = 1:size(tempBands,1)
        tempBands(ii,1) = bands(ii);
        tempBands(ii,2) = bands(ii+1) - dF;
    end
    bands = round(tempBands/dF);
    bands(end,2) = bands(end,2)+1;  % No -dF for the last one
else % error
    error('MATLAB:MEGanalysis:Improperparameter',...
        'bands must be either Nx1 or Nx2 array')
end
if bands(1,1)<1, bands(1,1)=1; end  % donot start from 0
% test if enough data in each band
while min(bands(:,2)-bands(:,1))<2*numRef
    b1 = bands(1,:);
    bands(2,1) = b1(1);
    bands(1,:)=[];
    if ~findNewCoefs  % band structure had ti be changed cannot use old 
                      % cleaning factors
        findNewCoefs=true;
        warning('MEGanalysis:notEnoughData', 'Cannot use old cleaning factors as bands had to be changed')
    end
    if size(bands,1)<3
        warning('MEGanalysis:minData','notEnough data for cleaning by bands')
    end
end
numBands = size(bands,1);

allCleaned=zeros(size(MEG));
% normalize the REFs
REF = REF - repmat(mean(REF,2), 1,size(REF,2));
rFactor = max(max(REF));
REF = REF/rFactor;
% The ref's may not be linearly independent project them on principal
% components, use all except the least significant 1
if size(REF,1)>6
    if isempty(eigVec)
        [V,D] = eig(REF*REF');
        eigVec = V;
    else
        V = eigVec;
    end
    REF =REF'*V;
    REF=REF(:,2:end)';
else
    eigVec = eye(size(REF,1));
end
numRef = size(REF,1);
% convert to frequency domain
if rot90
    REFfft = zeros(2*numRef, length(REF));
else
    REFfft = zeros(size(REF));
end
if rot90
    for ii = 1:numRef
        kk = 2*(ii-1) +1;
        REFfft(kk,:) = fft(REF(ii,:));
        REFfft(kk+1,:) =  REFfft(kk,:)*1i; % rotatae 90 degrees
    end
else
    for ii = 1:numRef
        REFfft(ii,:) = fft(REF(ii,:));
    end
end
if findNewCoefs
    if rot90
        allCoefs = zeros(numBands,2*numRef,numMeg);
    else
        allCoefs = zeros(numBands,numRef,numMeg);
    end
end

%% clean all channels
for ii=1:numMeg
    x=MEG(ii,:);
   if findNewCoefs
        if sum(abs(x)>overFlowTh)>0 && findNewCoefs
            % huge artifacts
            warning ('MATLAB:MEGanalysis:improperData',...
                'Big artifacts for channel %d. NOT cleaned',ii);
            allCleaned(ii,:) = x;
        else
            [cleanA1, coefA1] = clean1ByFft(x/rFactor, REFfft, samplingRate, bands*dF);
        end
    else  % use the supplied coefs
        oldCoefA = allCoefs(:,:,ii);
        [cleanA1, coefA1] = clean1ByFft(x/rFactor, REFfft, samplingRate, bands*dF,oldCoefA);
    end
    allCleaned(ii,:) = cleanA1*rFactor;
    allCoefs(:,:,ii) = coefA1;
end

bands = dF*bands;
bands(:,2) = bands(:,2)+dF;
bands = round(bands);
coefStruct = struct('allCoefs',allCoefs ,'eigVec',eigVec ,'bands',bands);

return
