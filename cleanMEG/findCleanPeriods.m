function cleanPeriods = findCleanPeriods(fileName, noiseLimit, tStep, lineFreq)
% find periods with no hi-freq. noise
%  cleanPeriods = findCleanPeriods(fileName, noiseLimit, tStep, lineFreq);
%
% fileName   - the pdf data file name
% noiseLimit - peaks at high freq / MAD(lowfrequencies) above which the
%              section is considered noisy [default 180]
% tStep      - analysis progresse in tStep in sec. [default 2]
% lineFreq   - the line frequency in Hz. [default 50].  the harmonics of
%              lf and their vicinity is ignored.
%
% cleanPeriods - cell array with list of clean periods.  One cell per MEG
%                channel.  in each cell (1,:) is the start time of the
%                clean periods and (2,:) is the end of the clean period in
%                seconds.
%
% analysis is done on 10 sec pieces which are advanced in steps of tStep.
%
% Example:
%   >> fileName = 'I:\Data\MEG\data\schizo\c,rfhp0.1Hz';
%   >> cleanPeriods = findCleanPeriods(fileName);


%Dec-2011  MA

%% initialize
if ~exist('tStep', 'var'), tStep=[]; end
if isempty(tStep), tStep = 2; end
if ~exist('lineFreq', 'var'), lineFreq=50; end
if isempty(lineFreq), lineFreq=50; end
if ~exist('noiseLimit', 'var'), noiseLimit=[]; end
if isempty(noiseLimit), noiseLimit=180; end
pdf = pdf4D(fileName);
if ispc % adjust this according to the max memory in the system
    samplesPerPiece = 250000;
else
    samplesPerPiece = 500000;
end
hdr = get(pdf,'header');
samplingRate   = double(get(pdf,'dr'));
numEpochs = length(hdr.epoch_data);
% lastSample = double(hdr.epoch_data{1}.pts_in_epoch);
if numEpochs>1
    error ('MATLAB:MEGanalysis:ImproperData',...
        'For now createCleanFileE doesnot work with Epoched data, try createCleanFile.')
else
    lastSample = double(hdr.epoch_data{1}.pts_in_epoch);
    epoched = false;
%     epochs = [1,lastSample];
end
% chit = channel_index(pIn,'TRIGGER');
chi = channel_index(pdf,'meg');
chn = channel_name(pdf,chi);
[~, chiSorted] = sortMEGnames(chn,chi);
% chix = channel_index(pIn,'EXTERNAL');
% chirf = channel_index(pIn,'ref');
% chie = channel_index(pIn,'EEG');
numChans = length(chi);


%% divide total time to pieces
if samplesPerPiece>lastSample
    samplesPerPiece= lastSample;
    aPieceOnly = true;
else
    aPieceOnly = false;
end
if ~epoched
    if  ~aPieceOnly
        numPieces = ceil(lastSample/samplesPerPiece);
        samplesPerPiece= floor(lastSample/numPieces);
        startApiece = 1:samplesPerPiece:lastSample;
        stopApiece  = startApiece+samplesPerPiece;
        deltaEnd = lastSample-stopApiece(end);
        if deltaEnd<0
            stopApiece(end) = lastSample;
            last = stopApiece(end)-startApiece(end);
            if last < samplesPerPiece/4  %unite with the previous piece
                startApiece(end)=[];
                stopApiece(end-1) = [];
            end
        else
            error ('wrong division of data')
        end
    else
        numSamples =lastSample; 
        numPieces = 1;
        firstS = 1;
        startApiece = firstS:samplesPerPiece:lastSample;
        stopApiece  = startApiece+samplesPerPiece;
        deltaEnd = lastSample-stopApiece(end);
        if deltaEnd<=0
            stopApiece(end) = lastSample;
        else
            error ('wrong division of data')
        end
    end
end
numPieces = length(startApiece);
cleanPeriods = cell(1,numChans);

%% search clean periods in each piece
for ii = 1:numPieces
    T0 = startApiece(ii)/samplingRate;
    MEG = read_data_block(pdf, [startApiece(ii) stopApiece(ii)], chiSorted);
    OKperiods = findAllOKperiods(MEG, samplingRate, noiseLimit, tStep, lineFreq);
    for jj = 1:numChans
        cleanPeriods{jj} = [cleanPeriods{jj} T0+OKperiods{jj}];
    end
end
    
%% unify the pieces
for ii = 1:numChans
    strt = cleanPeriods{ii}(1,:);
    ends = cleanPeriods{ii}(2,:);
    I = find(strt(2:end)<=ends(1:end-1));
    if ~isempty(I)
        ends(I) = [];
        strt(I+1) = [];
    end
    cleanPeriods{ii}= [strt ; ends];
end


return

