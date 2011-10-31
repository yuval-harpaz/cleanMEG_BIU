function [cleanA, coefA, eigVec] = clean1ByFft(y, ref, samplingRate, bands, oldCoef)
% clean MEG channel be fft
% [cleanA, coefA, rFactor] = clean1ByFft(y, ref, samplingRate, bands);
%
% y            - data og one MEG channel see below for scaling:
% ref          - the reference sensors or their DFT.  If ref is the time
%                domain data they should be as is.  If ref is the fft 
%                of ref, y should be scaled up by rFactor before calling
%                and cleanA should then be scaled back down
% samplingRate - in samples per sec.
% bands        - frequency bands to work on [default 0-4,4-8,8-16,...]
%              WARNING the default is based on sampling rate of ~2000/s
% oldCoef      - the coeficients for the channel at each band in the frequency domain
% cleanA       - the cleaned MEG channel (in the time domain)
% coefA        - the coeficients for each band in the frequency domain
% rFactor      - scaling factor for the REF channels (when in time-domain
% eigVec       - eigen vectors of ref in columns (Only if ref is in time
%                domain
% NOTES
%  It is assumed that REF (or their DFT) and MEG are scaled up by the same
%  factor


% Nov 2008  MA
% UPDATES
%  DEC-2008 If no improvement - leave the original data.  MA

%% initialize
if ~exist('oldCoef','var'), oldCoef=[]; end
if isempty(oldCoef)
    findNewCoef = true;
else
    coefA = oldCoef;
    findNewCoef = false;
end

y = y-mean(y);
[numRef, numData] = size(ref);
% if numRef>22
%     warning('MEGanalysis:tooManyRef', 'Only first 22 REF channels Used')
%     ref(23:end,:)=[];
%     numRef=22;
% end
if any(any(imag(ref)~=0)), refIsFreq=true;
else refIsFreq=false; end
if ~refIsFreq
   ref = ref-repmat(mean(ref),numRef,1);
   rFactor = max(max(ref));
   ref = ref/rFactor;
   % The ref's are not linearly independent project them on principal
   % components, use all except the least significant 1
   V = eig(ref*ref');
   eigVec=V;
   ref =ref'*V;
   ref=ref(:,2:end)';
   numRef = size(ref,1);
else
    rFactor=1;
    eigVec = eye(numRef);
end
y = y/rFactor;
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
    % truncate hi ranges if overlaps
    while bands(end,1)>=bands(end,2)-10
        bands(end,:)=[];
    end
    bands(end,2) = round((0.9*samplingRate/2-dF)/dF);
elseif bands(2,1)==bands(1,2)   % dF has to be subtracted
    if bands(3,1)<1/dF  % needs to be converted to dF units
        bands = round(bands/dF);
    end
    bands(2,:) = bands(2,:)-1;  % we assume here that bands is not in cycles but in dF units
else
    bands = round(bands/dF);
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
numBands = size(bands,1);

Y=fft(y);
if refIsFreq
    REF = ref;
else
    REF = zeros(size(ref));
    for ii = 1:numRef
        REF(ii,:) = fft(ref(ii,:));
    end
end
% cleanA = y;
if findNewCoef
    coefA = zeros(numBands,numRef);
end
CLEANa = zeros(size(y));

%% clean the fft in bands
% warning ('off', 'MATLAB:nearlySingularMatrix')
for bNo=1:numBands
    indx = bands(bNo,1):bands(bNo,2);
    Yp = Y(indx);
    C = REF(:,indx);
    if findNewCoef
        % A=C*C';
        % B=eye(numRef)/A;
        x = (Yp*C')/(C*C'); %*B;
        coefA(bNo,:)=x;
    else
        x= coefA(bNo,:);
    end
    CLEANa(indx) = Yp-x*C;
    if sum(CLEANa(indx)*CLEANa(indx)')>=sum(Yp*Yp') % no change needed
        CLEANa(indx) = Yp;
    end
end
% warning ('on', 'MATLAB:nearlySingularMatrix')

%%  smooth the transitions between bands
if Iseven(numData)
    mid = numData/2+1;
    CLEANa(mid+1:end) = conj(fliplr(CLEANa(2:mid-1)));
else
    mid = round((numData+1)/2);
    CLEANa(mid+1:end) = conj(fliplr(CLEANa(2:mid)));
end

%% convert back to time domain
cleanA = real(ifft(CLEANa))*rFactor;
% if ~refIsFreq
%     eigVec = eigVec*rFactor;
% % else
% %     coefA =coefA*rFactor; % scale back
% end
% coefA = coefA;  %  ??? *yFactor;

return
