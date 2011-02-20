function selectChans = selectChansForHB(MEG, samplingRate, chans2ignore)
% select channels with not to high a variance
%     selectChans = selectChansForHB(MEG, samplingRate, chans2ignore);
%
% MEG          - matrix of (Nchannels x Msamples)
% samplingRate - in samples per sec.
% chans2ignore - list of channels to be ignored in the analysis
%
% selectChans - list of channels that can be used for finding the QRS
% complex of EKG

% Feb-2011  MA

%% initialize
numChans = size(MEG,1);
minTooNoisy = round(numChans/25);
maxTooNoisy = round(numChans/5);
if ~exist('chans2ignore','var'), chans2ignore= []; end
chans2use = 1:numChans;
chans2use(chans2ignore) = [];
mf = MEG;
%% filter the data to eliminate slow drifts and fast noise
BandPassSpecObj=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
2,4,50,100,60,1,60,samplingRate);
BandPassFilt=design(BandPassSpecObj  ,'butter');
for ii = chans2use
    mf(ii,:) = myFilt(MEG(ii,:),BandPassFilt);
end
mf(chans2ignore,:)=0;
if max(max(mf)) <1e-7
    mf = mf*1e9;  % scale up to avoid underflow
end

%% lint away big chans
vr = var(mf,[],2);
mn = median(vr);
md = mad(vr);
tooNoisy = find(vr>(mn + 3*md));

%% test that reasonable
while length(tooNoisy)<minTooNoisy  % too few
    mf(tooNoisy,:) = 0;
    vr = var(mf,[],2);
    mn = median(vr);
    md = mad(vr);
    tooNoisy = unique([tooNoisy; find(vr>(mn + 3*md))]);
end
if length(tooNoisy) > maxTooNoisy  % too many
    [sortedVr,I] = sort(vr, 'descend');
    tooNoisy = I(1:maxTooNoisy);
end

selectChans = 1:248;
selectChans(unique([tooNoisy; chans2ignore']))=[];

return


