function noStepData = removeStep(Data, where, samplingRate, binSize, sec2smooth)
% Remove a step-like artifact from Data
%  noStepData = removeStep(Data, where, samplingRate, binSize, sec2smooth);
%
% Data         - Vector with data
% where        - where is the peak immediately after the step (see
%                findBigStep)
% samplingRate - in samples per second
% binSize      - fraction of sampling rate used to construct a smoothing
%                bin
% sec2smooth   - How many seconds to use for smoothing
%
% noStepData   - Data after removing the step.

% Dec-2010  MA
% UPDATES
%  Feb-2011 flip flag added to call of myConv.  MA

%% initialize
if ~exist('binSize', 'var'), binSize =[]; end
if isempty(binSize),binSize  = 1/10; end
if ~exist('sec2smooth', 'var'), sec2smooth =[]; end
if isempty(sec2smooth), sec2smooth = 15; end
noStepData = Data;
numData = length(Data);
madD = mad(Data);
bin = gaussBin(ceil(samplingRate*binSize),3);
I1 = where+ceil(sec2smooth*samplingRate);
if I1>numData, I1=numData; end
piece2smooth = Data(where:I1);
minMeet = round(2*(sec2smooth*samplingRate)/3);

%% smooth the piece
smoothedPiece = myConv(piece2smooth, bin, [], 'normalize','flip',false);
% add a linear piece at the very beginning
whereMet = find(abs(piece2smooth-smoothedPiece)<madD/20,1);
if whereMet>minMeet
    P = polyfit(10:whereMet, piece2smooth(10:whereMet), 1);
    aproxLine = polyval(P,1:whereMet);
    smoothedPiece(1:whereMet) = aproxLine;
end

%% clean the step
noStepData(where:I1) = noStepData(where:I1) - smoothedPiece;
% reduce the "spike" at the transient
noStepData(where-3:where) = noStepData(where);

return
