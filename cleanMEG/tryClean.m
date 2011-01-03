function [doLineF, doXclean, doHB, figH, QRS] = tryClean(MEG, samplingRate, Trig, XTR, doLineF, doXclean, doHB, chans2ignore, stepDur)
% Try to clean a piece of MEG to see if all can work
%   [doLineF, doXclean, doHB] = tryClean(MEG, samplingRate, Trig, XTR, ...
%                               doLineF, doXclean, doHB);
% 
% MEG          - the MEG channels
% samplingRate - in samples per s.
% Trig         - the trig channel
% XTR          - the external channels
% stepDur      - how long at top of a step the smoothing does apply
%                [default 0.2], but it depends on HiPass of MEG
%     The variable sbelow ar eboth input and output variables:
% doLineF, doXclean, doHB - flags (true/false) saying which cleaning method
%                you wish to apply or which are possible to apply
% 
% figH         - the fig handle for the plot of mean HB
% QRS          - the shape of the mean filtered QRS

% doLineF, doXclean, doHB - returned true if the requested cleaning can be
%                done
% chans2ignore  - list of channels with problems
% NOTE
%      at the beginning of this procedure there are a list of default
%      parameters which you mey need to change according to the machine
%      from which the data was collected.

% Dec-2010  MA

%% initialize
% some default params
lf = 50;                  % expected line frequrncy
Adaptive = true;          % how to clean the LF
Global = false;           % how to clean the LF
PhasePrecession = false;  % how to clean the LF
chans2ignore =  [74,204]; % known bad channels
outLierMargin = 10;       % channel with values that exceede 10 std are bad
hugeVal  = 1e-9;          % impossibly large MEG data
minVarAcc = 1e-4;         % XTR must have this variance or more
maxF = ceil(0.8*samplingRate/2);
if maxF>=140
    xBand = [1:139,140:20:maxF,maxF]'; % the XTR cleaning bands
    if xBand(end-1)==xBand(end)
        xBand(end)=[]; % avoid duplicaates
    end
else
    xBand = 1:maxF;    % the XTR cleaning bands
end
xChannels = 4:6;        % channels on which acceleration is recorded
HBperiod = 0.8;           % expected heart beat period in s.

% define missing params
if ~exist('chans2ignore', 'var'), chans2ignore = []; end
% if isempty (chans2ignore), chans2ignore = []; end
if ~exist('doLineF', 'var'), doLineF = []; end
if isempty (doLineF), doLineF = false; end
if ~exist('doXclean', 'var'), doXclean = []; end
if isempty (doXclean), doXclean = false; end
if ~exist('doHB', 'var'), doHB = []; end
if isempty (doHB), doHB = false; end
if ~exist('Trig', 'var'), Trig = []; end
if isempty (Trig), doLineF = false; end
if ~exist('XTR', 'var'), XTR = []; end
if isempty (XTR), doXclean = false; end
if ~exist('stepDur', 'var'), stepDur = []; end
if isempty (stepDur), stepDur = 1/5; end

numMEGchans = size(MEG,1);
chans2analyze = 1:numMEGchans;
chans2analyze(chans2ignore)=[];

%% remove big steps if any
bigStep=false(1,numMEGchans);
whereBig = nan(1,numMEGchans);
for iii=chans2analyze
    x = MEG(iii,:);
    [where, Tilda, atEnd] = findBigStep(x,samplingRate,[], stepDur, []);
    if ~isempty(where) && ~atEnd
        bigStep(iii)=true;
        whereBig(iii) = where;
    end
end
numBigSteps = sum(bigStep);
if numBigSteps>3 % clean the steps
    chans2clean = find(bigStep);
    for iii=chans2clean
        x = MEG(iii,:);
        where = whereBig(iii);
        MEG(iii,:) = removeStep(x, where, samplingRate);
    end
elseif numBigSteps>0 % marks this channels to be excluded
    chans2ignore = unique([chans2ignore find(bigStep)]);
end

%% check if any channel is too big or too noisy
ignoreChans = false(1, numMEGchans);
ignoreChans(chans2ignore) = true;
for iii=chans2analyze % 1:size(MEG,1)
    if sum(abs(MEG(iii,:))>hugeVal)>100 % there are >100 huge values
        warning('MEGanalysis:overflow',...
            ['In channel #' num2str(iii) ' more then 100 oveflows - ignored'])
        ignoreChans(iii)=true;
    end
end
vv = var(MEG,[],2);
outVv = (vv-median(vv))/mad(vv);
mVv = max(outVv);
if mVv > outLierMargin
    ovIndx = find(outVv>outLierMargin);
    for jj = 1:length(ovIndx)
        iii = ovIndx(jj);
        if ~ignoreChans(iii)  % add to ignored list
            warning('MEGanalysis:tooBig',...
                ['In channel #' num2str(iii) ' Variance is too big - ignored'])
            ignoreChans(iii)=true;
        end
    end
end
% sometimes the last value is HUGE??
[junkChan,junkData] = find(MEG(~ignoreChans,:)>hugeVal,1);
if ~isempty(junkData)
    endI = startI + size(MEG,2)-1;
    warning('MATLAB:MEGanalysis:nonValidData', ...
        ['MEGanalysis:overflow','Some MEG values are huge at: ',...
        num2str(junkData(1)/samplingRate) ' - truncated'])
    MEG(junkChan,junkData:end)=0;
    ignoreChans(junkChan)=true;
end

chans2analyze(chans2ignore)= [];
MEG(chans2ignore,:)=0;

%% find the LF bit in Trig.  Try 1st the 256 bit
if doLineF
    lf1 = samplingRate/lf;
    bit = uint16(256);
    trig= uint16(Trig);
    f=bitand(trig,bit);
    whereUp = find(diff(f)==bit);
    meanLF =mean(diff(whereUp));
    if ~((lf1-2)<meanLF && meanLF<(lf1+2)) %  not within +/- 2 samples
        whichLF = findLFbit(trig,samplingRate);  % find the LF bit
        if isempty(lineF)
            warning('MATLAB:MEGanalysis:noData','Line Frequency trig not found')
            doLineF = false;
        else
            whereUp = find(diff(f)==whichLF);
        end
    end
end

%% clean the Line Frequency
if doLineF
    LFcycleLength = max(diff(whereUp)+3);
    % meanPeriod = mean(whereUp);
    if Global || Adaptive  % use one average for the entire file
        MEGlfCycle = zeros(numMEGchans,LFcycleLength);
        if exist('XTR','var')
            XTRlfCycle = zeros(size(XTR,1),LFcycleLength);
        end
    elseif PhasePrecession % use interpolation
        MEGlfCycle = zeros(numMEGchans,3*LFcycleLength);
        if exist('XTR','var')
            XTRlfCycle = zeros(length(chixSorted),3*LFcycleLength);
        end
    else
        disp(['The only legal cleaning methods are: ''Global'''...
            ', ''Adaptive'', or ''PhasePrecession'''])
    end
    totalCycles =0; % cycles of LF
    
    %% check if any channel is too big or too noisy
    ignoreChans = false(1, numMEGchans);
    ignoreChans(chans2ignore) = true;
    for iii=chans2analyze % 1:size(MEG,1)
        if sum(abs(MEG(iii,:))>hugeVal)>100 % there are >100 huge values
            warning('MEGanalysis:overflow',...
                ['In channel #' num2str(iii) ' more then 100 oveflows - ignored'])
            MEG(iii,:)=0;
            ignoreChans(iii)=true;
        end
    end
    vv = var(MEG,[],2);
    outVv = (vv-median(vv))/mad(vv);
    mVv = max(outVv);
    if mVv > outLierMargin
        ovIndx = find(outVv>outLierMargin);
        for jj = 1:length(ovIndx)
            iii = ovIndx(jj);
            if ~ignoreChans(iii)  % add to ignored list
                warning('MEGanalysis:tooBig',...
                    ['In channel #' num2str(iii) ' Variance is too big - ignored'])
                % MEG(iii,:)=0;
                ignoreChans(iii)=true;
            end
        end
    end
    % sometimes the last value is HUGE??
    [junkChan,junkData] = find(MEG(~ignoreChans,:)>hugeVal,1);
    if ~isempty(junkData)
        endI = startI + size(MEG,2)-1;
        warning('MATLAB:MEGanalysis:nonValidData', ...
            ['MEGanalysis:overflow','Some MEG values are huge at: ',...
            num2str(junkData/samplingRate) ' - truncated'])
        MEG(junkChan,junkData:end)=0;
        ignoreChans(junkChan)=true;
    end
    
    %% find the the position and mean of LF cycle
    % find where are the LF cycles
    numCyclesHere = length(whereUp);
    if Global % find the overall mean LF cycle
        maxL = LFcycleLength-1;
        sumC = zeros(numMEGchans,maxL+1);
        for cycle = 1:numCyclesHere-1;
            startCycle = whereUp(cycle);
            if startCycle+maxL <= size(MEG,2)
                sumC = sumC + MEG(:,startCycle:startCycle+maxL);
            end
        end
        MEGlfCycle = MEGlfCycle+sumC;
        totalCycles = totalCycles+numCyclesHere-1;
        sumC = zeros(numREFchans,maxL+1);
        for cycle = 1:numCyclesHere-1;
            startCycle = whereUp(cycle);
            if startCycle+maxL <= size(MEG,2)
                sumC = sumC + REF(:,startCycle:startCycle+maxL);
            end
        end
        if exist('XTR','var')
            sumC = zeros(size(XTR,1),maxL+1);
            for cycle = 1:numCyclesHere-1;
                startCycle = whereUp(cycle);
                if startCycle+maxL <= size(MEG,2)
                    sumC = sumC + XTR(:,startCycle:startCycle+maxL);
                end
            end
            XTRlfCycle = XTRlfCycle+sumC;
        end
    end % end of finding global LF cycles
    
    if Global
        MEGlfCycle = MEGlfCycle/totalCycles;
        MEGlfCycle = MEGlfCycle - repmat(mean(MEGlfCycle,2),1,maxL+1);
        if exist('XTR','var')
            XTRlfCycle = XTRlfCycle/totalCycles;
            XTRlfCycle = XTRlfCycle - repmat(mean(XTRlfCycle,2),1,maxL+1);
        end
    end % end of normalizing global averages
%     if Adaptive % oreoare arrays for ongoing mean LF cycles
%         meanL = max(diff(whereUp)) +1;
%         meanMEG = zeros (numMEGchans,meanL);
%         if exist('XTR','var')
%             meanXTR = zeros (size(XTR,1),meanL);
%         end
%     end
    
    %%    % clean the 50Hz if needed
    for nc = chans2analyze
        x=MEG(nc,:);
        if Global
            [y,Tilda] = cleanMean(x, whereUp, MEGlfCycle(nc,:), 1, []);
        elseif Adaptive
            [y, Tilda] = cleanLineF(x, whereUp, [], 'Adaptive');
            % define the array of means
%         elseif phasePrecession % must be phase precession
%             y = cleanLineF(x, whereUp, [], 'phasePrecession');
        else
            error ('MATLAB:MEGanalysis:IllegalParam'...
                ,['Allwed METHODs are: ' legalArgs])
        end
        MEG(nc,:)=y;
    end
    if exist('XTR','var') && ~isempty(XTR) % clean the 50 Hz
        for nc = 1:size(XTR,1)
            x=XTR(nc,:);
            if Global
                y = cleanMean(x, whereUp, XTRlfCycle(nc,:), 1, []);
            elseif Adaptive
                [y, Tilda] = cleanLineF(x, whereUp, [], 'Adaptive');
                % define the array of means
%             elseif phasePrecession % must be phase precession
%                 y = cleanLineF(x, whereUp, [], 'phasePrecession');
            else
                error ('MATLAB:MEGanalysis:IllegalParam'...
                    ,['Allwed METHODs are: ' legalArgs])
            end
            XTR(nc,:)=y;
        end
        % check that there is a signal on xChannels
        if size(XTR,1)<max(xChannels)  % requested XTR does not exist
            doXclean = false;
        elseif any(var(XTR(xChannels,:),[],2)<minVarAcc)
            doXclean=false;
        end
    end
end

%% clean the XTR channels
if doXclean
    for nc = chans2analyze
        x=MEG(nc,:);
        y = coefsAllByFFT(x, XTR(xChannels,:),...
            samplingRate, xBand, [], 1,false);
        MEG(nc,:) = y;
    end
end

%% clean the Heart Beat
if doHB
    mMEG = mean(MEG,1);
    [whereisHB, zTime, Errors, amplitudes, meanBeat, figH, QRS]= findHB01(mMEG, samplingRate,HBperiod,...
        'PLOT', 'VERIFY');
    drawnow
    if (sum(Errors.shortHB)>3) || (sum(Errors.longHB)>3) || (Errors.numSmallPeaks>3) % added by Yuval
        warning('MATLAB:MEGanalysis:ImproperData',...
            'Heart beat signal not clear enough No HB cleaning done');
        disp(Errors)
        doHB = false;
        return
    end % end of added by Yuval
    % code from here on are pieces used to actually clean the HB - but it
    % was not debugged, and is not needed for the purpose of the tryClean
    % procedure
%     HBcycleLength = round(HBperiod*samplingRate);
%     % from here on zTime, iBefore and iAfter are FIXED!
%     iBefore = zTime;
%     iAfter = HBcycleLength -iBefore-1; 
%     amplitudes = amplitudes/mean(amplitudes);
%     if sum(amplitudes<=0)  > 0
%         doHB= false;
%         return
%     end
%     mMEGhb = zeros(numMEGchans,HBcycleLength);
%     for chan = 1:numMEGchans
%         x=MEG(chan,:);
%         mMEGhb(chan,:) = meanAround(x, whereisHB,iBefore, iAfter);
%     end
%     cycleLength = iBefore+iAfter+1;
%     % offset so edges are near zero
%     for jj = 1:size(mMEGhb,1)
%         tailOffset = mean([mMEGhb(jj,1:10) mMEGhb(jj,end-9:end)]);
%         mMEGhb(jj,:) = mMEGhb(jj,:) - tailOffset;
%     end
%     numChans = size(mMEGhb,1);
%     numBeats = length(whereisHB);
%     % cycleLength = size(mMEGhb,2);
%     % Use only the large QRS shape to normalize amplitude.
%     nTmplt = MEGhbCycle(:, QRSstartI:QRSendI);
%     nTmplt = nTmplt-repmat(mean(nTmplt,2),1,size(nTmplt,2));
%     for jj=1:numChans
%         x= nTmplt(jj,:);
%         nTmplt(jj,:) = x/sqrt(x*x');
%     end
%     nTmplt = nTmplt'; % each normalized template in one column
%     
%     % test if first/lst not too close to th edge
%     if whereisHB(1)-QRSbefore <=0
%         jStrt=2;
%     else
%         jStrt=1;
%     end
%     if whereisHB(end)+QRSafter >size(MEG,2)
%         jEnd=numBeats-1;
%     else
%         jEnd=numBeats;
%     end
%     
%     for chan = 1:numChans
%         Amplitudes = zeros(size(whereisHB));
%         x=MEG(chan,:);
%         for jj = jStrt:jEnd  %BUT take care of first and last cycles!!
%             % what if too near the begining? =XXX
%             thisBeat = x((whereisHB(jj)-QRSbefore):(whereisHB(jj)+QRSafter));
%             thisBeat = thisBeat-mean(thisBeat);
%             Amplitudes(jj) = thisBeat*nTmplt(:,chan);
%         end
%         Amplitudes(jStrt:jEnd) = Amplitudes(jStrt:jEnd)/mean(Amplitudes(jStrt:jEnd));
%         if jStrt>1
%             Amplitudes(1)=1;
%         end
%         if jEnd>numBeats
%             Amplitudes(end)=1;
%         end
%         % clip amplitudes
%         if clipAmplitude
%             Amplitudes(Amplitudes<minAmplitude) = minAmplitude;
%             Amplitudes(Amplitudes>maxAmplitude) = maxAmplitude;
%         end
%         allAmplMEG(chan,1:length(Amplitudes)) = Amplitudes;
%         y = cleanMean(x, whereisHB, MEGhbCycle(chan,:), zTime, Amplitudes);
%         MEG(chan,:)=y;
%     end
else
figH='no heart fig';QRS='no qrs';
end

return

