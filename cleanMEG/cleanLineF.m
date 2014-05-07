function [cleaned, mean1,noiseSamp] = cleanLineF(dataA, whereUp, epochs, method, mean0,startNum,noiseThr)
%  clean the line frequency based on points at which the Mains flipped from
%  negative to positive
%    [cleaned, mean1] = cleanLineF(dataA, whereUp, epochs, method, Mean0);
% dataA   - array with one channel of MEGdata
% whereUp - list of column indices where Line Freq changed from negative to
%           positive.  If the mark of this change comes from ExtTrig #1 it
%           may be found by:
%                           whereUp=find(diff(mod(trig5to45,511)>=256)==1);
% epochs  - list of indices where the last piece ended [default [0,length(dataA)]
% method  - one of 3 possibilities:
%            'GLOBAL' - the average over the entire data is used [default]
%          'ADAPTIVE' - The averge gradually follow the mean, starting with
%                       the average over 5 sec.
%   'PHASEPRECESSION' - The data is interpolated 10 folds and the phase
%                       precession at each trig is considered
% Mean0   - start from that mean.  If not available comput from 1-st 256
%          cycles.
% startNum- how many cycles to take for adaptive template
%
% cleaned - same dimension like dataA but cleaned from the line frequency
%           artefact.
% mean1   - the last mean used
%
% noiseThr- how many std(mean(abs(data of one cycle))) to consider as
%           noise. default is 5 SD.

%  Sep-2008  MA
% UPDATES
%  Dec-2008  cycle end is considered by meanL and not maxL  MA
%             Adapted for data composed of non-contiguos pieces
%  Oct-2010  3 methods of cleaning added
%  Feb-2011  Bugs in adaptive methods fixed
%  Apr-2014  leave low frequencies in the data, reject noisy segments from
%            template, allow choose startNum and more. see
%            github/yuval-harpaz for all changes

%% initialize
if nargin>3
    okArgs = {'GLOBAL','ADAPTIVE','PHASEPRECESSION'};
    k = find(strncmpi(method, okArgs,6));
    if isempty(k)
        error('MATLAB:MEGanalysis:BadParameter',...
            'Unknown method name:  %s.',method);
    elseif length(k)>1
        error('MATLAB:MEGanalysis:BadParameter',...
            'Ambiguous method name:  %s.',method);
    else
        switch(k)
            case 1  % GLOBAL
                Global=true;
                Adaptive = false;
                phasePrecession = false;
            case 2 % ADAPTIVE
                Global=false;
                Adaptive = true;
                phasePrecession = false;
                if ~exist('startNum','var')
                    startNum=[];
                end
                if isempty(startNum)
                    startNum=256;
                end
            case 3 % PhasePrecession
                Global=false;
                Adaptive = false;
                phasePrecession = true;
        end
    end
else
    Global=true; % the default
    Adaptive = false;
    phasePrecession = false;
end
if ~exist('mean0', 'var'), mean0=[]; end
if ~isempty(mean0)
    maxL = length(mean0)-1;
else
    maxL = max(diff(whereUp));
end
if ~exist('epochs', 'var'), epochs = []; end
if isempty(epochs), epochs = [0, length(dataA)]; end
if length(epochs)>2
    epoched = true;
    numEpochs = length(epochs)-1;
    epochS = epochs(1:end-1)+1;
    epochE = [(epochS(2:end)-1) length(dataA)];
else
    epoched = false;
    numEpochs = 1;
    epochS = epochs(1)+1;
    epochE = epochs(2)-1;
    if length(whereUp)<500 && ~Global
        warning('MATLAB:MEGanalysis:ImproperCombination',...
            'Data is too short so only GLOBAL is allowed')
        Global=true; % the default
        Adaptive = false;
        phasePrecession = false;
    end
end
if epoched && ~Global
    warning('MATLAB:MEGanalysis:ImproperCombination',...
        'Data is epoched and only GLOBAL is allowed')
    Global=true; % the default
    Adaptive = false;
    phasePrecession = false;
end
% adjust if first epoch does not start at 1
if epochS(1)>1
    epochS=[1 epochS];
    epochE = [epochS(2)-1 epochE];
    numEpochs = length(epochS);
end
dW = diff(whereUp);
if ~exist('noiseThr','var')
    noiseThr=[];
end
if isempty(noiseThr)
    noiseThr=5;
end


%% Clean
if Global
    if ~epoched
        meanL=round(mean(diff(whereUp)));
        firstCycleStart = whereUp(1);
        first=1;
        last = find(whereUp<=epochE,1, 'last')-1;
        kk=1;
        meanLapprox = round(mean(diff(whereUp)));  % temporary evaluation of a cycle
        shortEpoch = false;
        noUpInEpoch = false;
        % lastCycleStart = whereUp(end);
    else % consider only whereUp with complete cycles at each epoch
        sumPeriods = 0;
        numPeriods = 0;
        firstCycleStart = ones(1,numEpochs);
        lastCycleStart  = ones(1,numEpochs);
        meanLapprox = round(mean(diff(whereUp)));  % temporary evaluation of a cycle
        %  check if last epoch contains less then 1 cycle
        vLast = find(whereUp>epochS(end),1);
        if isempty(vLast)
            numEpochs = numEpochs-1;
        end
        first = ones(1,numEpochs);
        last = ones(1,numEpochs);
        % find where is the first epoch after the first whereUp
        if whereUp(2)>epochE(1)
            shortEpoch=true;
            first(1) =1;
            last(1) = 1;
            kk=2;
            if whereUp(1)>epochE(1)
                noUpInEpoch=true;
                % find where in a cycle is the first epoch
            else
                noUpInEpoch=false;
            end
        else
            shortEpoch=false;
            kk=1;
        end
        for ii=kk:numEpochs
            first(ii) = find(whereUp>epochS(ii),1);
            ll = find(whereUp((first(ii)+1):end)<=epochE(ii),1, 'last');
            if ~isempty(ll)
                if isempty(last(ii)), last(ii)=length(whereUp);
                else
                    ll = ll+first(ii);
                    last(ii) = ll;
                end
                firstCycleStart(ii) = whereUp(first(ii));
                lastCycleStart(ii) = whereUp(last(ii));
                if last(ii)>length(dW), last(ii)=last(ii)-1; end
                sumPeriods = sumPeriods +sum(dW(first(ii):last(ii)));
                numPeriods =numPeriods + last(ii)-first(ii)+1;
            end
        end
        meanL = round(sumPeriods/numPeriods);
    end
    % lastCycleStart  = whereUp(end);
    cleaned = zeros(size(dataA));
    
    %% get the mean signal per line cycle
    meanLine = zeros(size(dataA,1),max(dW)+1);
    numData = zeros(1,max(dW)+1);
    lastDataSample = epochS(1)-1;
    noiseSamp=[];
    for ii=kk:numEpochs
        if epochE(ii)-epochS(ii)>3*meanLapprox  % do not use too short pieces for meanLine
            [mL,nSamp] = oneLineCycle(dataA(:,(epochS(ii):epochE(ii))), ...
                whereUp((first(ii):last(ii)))-lastDataSample);
            numInThisMean  = size(mL,2);
            numData(1:numInThisMean) = numData(:,1:numInThisMean) +...
                (last(ii) -first(ii) +1);
            meanLine(:,1:numInThisMean) = meanLine(:,1:numInThisMean) +...
                (last(ii) -first(ii) +1)*mL;
            lastDataSample = epochE(ii);
            noiseSamp=[noiseSamp,nSamp+epochS(ii)-1];
        end
    end
    extraI = find(numData==0,1);
    meanLine(:,extraI:end)=[];
    numData(extraI:end)=[];
    meanLine = meanLine./repmat(numData,size(meanLine,1),1);
    
    %% subtract from signal one cycle at a time
    % subtract the initial piece
    %NOTE if jj==1 and whereup strats just there there may be a problem
    if shortEpoch % treat the first one separately
        if noUpInEpoch % treat by the up in next
            nStart = whereUp(1)-epochS(1)+1;
            nEnd   = whereUp(1)-epochE(1) +1;
            cleaned(:,1:epochE(1)) = dataA(:,1:epochE(1)) - ...
                meanLine(:,(end-nStart):(end-nEnd));
        else % relate to the first up
            boundery = whereUp(1);
            nBefore = boundery-1;
            nAfter = epochE(1)-boundery +1;
            cleaned(:,1:nBefore) = dataA(:,1:nBefore)- ...
                meanLine(:,end-nBefore+1:end);
            cleaned(:,epochE(1)-nAfter+1:epochE(1)) = dataA(:,epochE(1)-nAfter+1:epochE(1))- ...
                meanLine(:,1:nAfter);
        end
    end
    for jj=kk:numEpochs
        numInThisCycle = firstCycleStart(jj) -epochS(jj);
        if meanL==numInThisCycle
            strtOffset=1;
        else
            strtOffset=0;
        end
        if epochE(jj)-whereUp(last(jj))>meanL  %  add one cycle
            thisWhereUp = whereUp(first(jj)-strtOffset:(last(jj)+1));
        else
            thisWhereUp = whereUp(first(jj)-strtOffset:last(jj));
        end
        if meanL>numInThisCycle
            cleaned(:,epochS(jj):firstCycleStart(jj)-1) = dataA(:,epochS(jj):firstCycleStart(jj)-1)...
                - meanLine(:,(meanL-numInThisCycle):meanL-1);
            iEnds = thisWhereUp(jj); %?????????????
        end
        % subtract from the center part
        for ii=1:length(thisWhereUp)-1
            iStrt = thisWhereUp(ii);
            iEnds = thisWhereUp(ii+1) -1;
            numInThisCycle = iEnds-iStrt+1;
            cleaned(:,iStrt:iEnds) = dataA(:,iStrt:iEnds)-meanLine(:,1:numInThisCycle);
        end
        % subtract the leftover tail
        if epochE(jj)>iEnds  ; %clean the tail
            lastHere = epochE(jj);
            cleaned(:,iEnds+1:lastHere) = dataA(:,iEnds+1:lastHere)-...
                meanLine(:,1:lastHere-iEnds);
            
        end
    end
    if lastHere<length(dataA)  % clean the last piece
        lastTail = size(cleaned,2)-lastHere;
        if lastHere-iEnds+lastTail<=meanL % less then one cycle left
            cleaned(:,lastHere+1:end) = dataA(:,lastHere+1:end)-...
                meanLine(:,lastHere-iEnds+1:lastHere-iEnds+lastTail);
        else % treat later
            warning('MATLAB:MEGanalysis:incompleteCalculations', ['Last ' num2str(lastTail) ' Not cleaned!']);
            cleaned(:,lastHere+1:end) = dataA(:,lastHere+1:end);
        end
    end
    mean1 = meanLine;
elseif Adaptive
    %% generate a slowly changing average
    cleaned = dataA;
    %    startNum=256;
    numCycles = length(whereUp);
    Q = 1-1/startNum;
    sum1 = zeros(1,maxL+1);
    ml1 = nan(numCycles,maxL +1);
    if ~exist('mean0', 'var')
        mean0 = [];
    else
        if sum(abs(mean0))==0, mean0=[]; end
    end
    % Estimate Noise
    for cycle = 1:(numCycles-2)
        startCycle = whereUp(cycle);
        amp1(cycle) = mean(abs(dataA(startCycle:startCycle+maxL)-mean(dataA(startCycle:startCycle+maxL)))); %#ok<AGROW>
%         if ~isempty(find([startCycle:startCycle+maxL]==688726))
%             display('noise')
%         end
    end
    amp2=(amp1-mean(amp1))./std(amp1);
    %cyci=find(amp2<noiseThr);
    noise=min(amp1(amp2>=noiseThr));
    if isempty(noise)
        noise=max(amp1); % to accept all segments
    end
    %% compute a simple average
    if isempty(mean0)  % compute for the first 256 (or startNum)
        %         for cycle = 1:startNum
        %             startCycle = whereUp(cycle);
        %             amp1(cycle) = mean(abs(dataA(startCycle:startCycle+maxL)-mean(dataA(startCycle:startCycle+maxL)))); %#ok<AGROW>
        %         end
        %         amp2=(amp1-mean(amp1))./std(amp1);
        %         cyci=find(amp2<noiseThr);
        %         noise=min(amp1(amp2>=noiseThr));
        %         if isempty(noise)
        %             noise=max(amp1); % to accept all segments
        %         end
        cycCount=0;
        noiseSamp=[];
        for cycle = 1:startNum
            startCycle = whereUp(cycle);
            if mean(abs(dataA(startCycle:startCycle+maxL)-mean(dataA(startCycle:startCycle+maxL))))<=noise
                sum1 = sum1 + dataA(startCycle:startCycle+maxL);
                cycCount=cycCount+1;
            else
                noiseSamp=[noiseSamp,startCycle:(startCycle+maxL)]; %#ok<*AGROW>
            end
        end
        ml1(1:startNum,:) = repmat(sum1/cycCount,startNum,1);
    else % mean0 was provided
        % check that OK
        r = size(mean0,1);
        if r==1 % a row vector - transpose
            mean0 = mean0';
            r = size(mean0,1);
        end
        if r~= length(sum1)
            error('MATLAB:MEGanalysis:ImproperParam',...
                'The initial mean must be %d long', length(sum1))
        end
        ml1(1:startNum,:) = repmat(mean0',startNum,1);
    end  % end of getting the first startNum averages
    % continue in adaptive way
    for cycle = startNum+1:numCycles
        startCycle = whereUp(cycle);
        if startCycle+maxL <= size(dataA,2)
            if mean(abs(dataA(startCycle:startCycle+maxL)...
                    -mean(dataA(startCycle:startCycle+maxL)))) <= noise
                ml1(cycle,:) = Q*ml1(cycle-1,:) + ...
                    dataA(startCycle:startCycle+maxL)/startNum;
            else
                noiseSamp=[noiseSamp,startCycle:(startCycle+maxL)];
                ml1(cycle,:)=ml1(cycle-1,:);
            end
        else % extra cycles copy the previous one
            ml1(cycle,:)=ml1(cycle-1,:);  % copy the last one
        end
    end
    % BL correction for template
    for tempi=1:size(ml1,1)
        ml1(tempi,:)=ml1(tempi,:)-mean(ml1(tempi,:));
    end
    for ii=1:length(whereUp)-1
        iStrt = whereUp(ii);
        iEnds = whereUp(ii+1) -1;
        numInThisCycle = iEnds-iStrt+1;
        cleaned(iStrt:iEnds) = dataA(iStrt:iEnds)-ml1(ii,1:numInThisCycle);
    end
    % treat the edges
    if whereUp(1)>1  %header before first whereUp
        numInHeader = whereUp(1)-1;
        cleaned(1:numInHeader) = dataA(1:numInHeader)...
            - ml1(1,end-numInHeader+1:end);
    end
    if whereUp(end)<length(dataA) % tail after whereUp
        numInTail = length(dataA)-whereUp(end);
        cleaned(end-numInTail:end) = dataA(end-numInTail:end)...
            - ml1(end, end-numInTail:end);
    end
    mean1 = ml1(end,:);
elseif phasePrecession
    interpNo =10; % How many interpolation points between samples
    [cleaned, mean1] = cleanWphase(dataA,whereUp,interpNo);
else
    error('MATLAB:MEGanalysis:unknownParam','method was not defined')
end
return