function [doLineF, doXclean, doHB, figH, QRS] = tryClean(MEG, samplingRate,...
    Trig, XTR, xChannels, doLineF, doXclean, doHB, chans2ignore, stepDur, hugeVal,ECG,HBperiod)
% Try to clean a piece of MEG to see if all can work
%    [doLineF, doXclean, doHB, figH, QRS, chans2ignore] = tryClean...
%           (MEG, samplingRate, Trig, XTR, doLineF, doXclean, doHB,...
%            stepDur, hugeVal);
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
% hugeVal      - What is considered as overflow
%
% figH         - the fig handle for the plot of mean HB
% QRS          - the shape of the mean filtered QRS
% doLineF, doXclean, doHB - returned true if the requested cleaning can be
%                done
% chans2ignore  - list of channels with problems
% NOTE
%      at the beginning of this procedure there are a list of default
%      parameters which you may need to change according to the machine
%      from which the data was collected.

% Dec-2010  MA

%% initialize
% some default params
lf = 50;                  % expected line frequrncy
Adaptive = lf==lf;        % how to clean the LF
Global = ~Adaptive;       % how to clean the LF
outLierMargin = 20;       % channel with values that exceede 20 std are bad
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
% xChannels = 4:6;        % channels on which acceleration is recorded
% HBperiod = 0.8;           % expected heart beat period in s.

% define missing params
if ~exist('chans2ignore', 'var'), chans2ignore = []; end
if isempty (chans2ignore), chans2ignore = [74, 204]; end
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
if ~exist('hugeVal','var'), hugeVal=[]; end
if isempty(hugeVal)
    hugeVal  = 1e-9;      % impossibly large MEG data
end

figH=[];  % in case not HB is done
QRS=[];

numMEGchans = size(MEG,1);
ignoreChans = false(1,numMEGchans);
ignoreChans(chans2ignore)=true;
chans2analyze = find(~ignoreChans);

%% remove big steps if any
bigStep=false(1,numMEGchans);
whereBig = nan(1,numMEGchans);
for iii=chans2analyze
    x = MEG(iii,:);
    [where, ~, atEnd] = findBigStep(x,samplingRate,[], stepDur, []);
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
elseif numBigSteps>0 % marks this channels to be excluded for now
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
    warning('MATLAB:MEGanalysis:nonValidData', ...
        ['MEGanalysis:overflow','Some MEG values are huge at: ',...
        num2str(junkData(1)/samplingRate) ' - truncated'])
    MEG(junkChan,junkData:end)=0;
    ignoreChans(junkChan)=true;
end

chans2analyze=find(~ignoreChans);
MEG(ignoreChans,:)=0;

%% find the LF bit in Trig.  Try 1st the 256 bit
if doLineF
    bit = findLFbit(Trig,samplingRate);
    if ~isempty(bit)
        trig= uint16(Trig);
        f=int16(bitand(trig,bit));
        whichLF = int16(bit);
        whereUp = find(diff(f)==whichLF);
    else
        warning('MATLAB:MEGanalysis:noData','Line Frequency trig not found')
        doLineF = false;
    end
end

%% clean the Line Frequency
if doLineF
    LFcycleLength = max(diff(whereUp)+3);
    % meanPeriod = mean(whereUp);
    if Global || Adaptive  % use one average for the entire file
        if exist('XTR','var')
            XTRlfCycle = zeros(size(XTR,1),LFcycleLength);
        end
    end
    
    %% check if any channel is too big or too noisy
    ignoreChans = false(1, numMEGchans);
    ignoreChans(chans2ignore) = true;
    chans2analyze = find(~ignoreChans);
%     for iii=chans2analyze % 1:size(MEG,1)
%         if sum(abs(MEG(iii,:))>hugeVal)>100 % there are >100 huge values
%             warning('MEGanalysis:overflow',...
%                 ['In channel #' num2str(iii) ' more then 100 oveflows - ignored'])
%             MEG(iii,:)=0;
%             ignoreChans(iii)=true;
%         end
%     end
%     vv = var(MEG,[],2);
%     outVv = (vv-median(vv))/mad(vv);
%     mVv = max(outVv);
%     if mVv > outLierMargin
%         ovIndx = find(outVv>outLierMargin);
%         for jj = 1:length(ovIndx)
%             iii = ovIndx(jj);
%             if ~ignoreChans(iii)  % add to ignored list
%                 warning('MEGanalysis:tooBig',...
%                     ['In channel #' num2str(iii) ' Variance is too big - ignored'])
%                 % MEG(iii,:)=0;
%                 ignoreChans(iii)=true;
%             end
%         end
%     end
    % sometimes the last value is HUGE??
    [junkChan,junkData] = find(MEG(~ignoreChans,:)>hugeVal,1);
    if ~isempty(junkData)
%         endI = startI + size(MEG,2)-1;
        warning('MATLAB:MEGanalysis:nonValidData', ...
            ['MEGanalysis:overflow','Some MEG values are huge at: ',...
            num2str(junkData/samplingRate) ' - truncated'])
        MEG(junkChan,junkData:end)=0;
        ignoreChans(junkChan)=true;
%         chans2ignore = find(ignorechans);
        chans2analyze = find(~ignoreChans);
    end
    
    %% find the the position and mean of LF cycle
    % find where are the LF cycles
    numCyclesHere = length(whereUp);
    maxL = LFcycleLength-1;
    MEGlfCycle = zeros(length(chans2analyze), maxL);
    totalCycles = 0;
    if Global % find the overall mean LF cycle
        sumC = zeros(numMEGchans,maxL+1);
        for cycle = 1:numCyclesHere-1;
            startCycle = whereUp(cycle);
            if startCycle+maxL <= size(MEG,2)
                sumC = sumC + MEG(:,startCycle:startCycle+maxL);
            end
        end
        MEGlfCycle = MEGlfCycle+sumC;
        totalCycles = totalCycles+numCyclesHere-1;
        if exist('XTR','var')
            sumC = zeros(size(XTR,1),maxL+1);
            XTRlfCycle = zeros(size(sumC));
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
    
    %%    % clean the 50Hz if needed
    for nc = chans2analyze
        x=MEG(nc,:);
        if Global
            [y,~] = cleanMean(x, whereUp, MEGlfCycle(nc,:), 1, []);
        elseif Adaptive
            [y, ~] = cleanLineF(x, whereUp, [], 'Adaptive');
        end
        MEG(nc,:)=y;
    end
    if exist('XTR','var') && ~isempty(XTR) % clean the 50 Hz
        for nc = 1:size(XTR,1)
            x=XTR(nc,:);
            if Global
                y = cleanMean(x, whereUp, XTRlfCycle(nc,:), 1, []);
            elseif Adaptive
                [y, ~] = cleanLineF(x, whereUp, [], 'Adaptive');
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
    
    mMEG = mean(MEG(chans2analyze,:),1);
    if exist('ECG','var');
        if ~isequal(size(ECG),size(mMEG))
            error('ECG channel not the right size')
        else
            mMEG=ECG;
        end
    end
        
        
    [~, ~, Errors, ~, ~, figH, QRS]= findHB01(mMEG, samplingRate,HBperiod,...
        'PLOT', 'VERIFY');
    drawnow
    if (length(Errors.shortHB)>3) || (length(Errors.longHB)>3) || (Errors.numSmallPeaks>3) % added by Yuval
        warning('MATLAB:MEGanalysis:ImproperData',...
            'Heart beat signal not clear enough No HB cleaning done');
        disp(Errors)
        doHB = false;
        return
    end % end of added by Yuval
end

% chans2ignore = find(ignoreChans);
return

