function OKperiods = findAllOKperiods(MEG, samplingRate, noiseLimit, tStep, lineFreq)
%  Find 10 sec periods without big spectral lines at hi frequency
%    noHiFnoise = findOKperiods(x, samplingRate, noiseLimit, tStep,...
%    lineFreq);
%
% MEG          - the MEG data [Ncahnnels X Mdatapoints]
% samplingRate - in samples per second
% noiseLimit   - max allowed ratio of hi-f peaks relative to mean of lower
%                frequencies [default 100]
% tStep        - steps to advance the analysis, in seconds [default 5 sec]
% lineFreq     - [default 50]
%
% OKperiods - cell array of 1xNchannels matrices, each of 2xN elements
%             with start and end points of periods with no high
%             frequency noise.
%
%   ALGORITHM
% The PSD of pieces 10 s long is computed.  The median and mad of
% frequencies between 103 and 147 Hz is computed.  If the highest peak for
% frequencies above 160 Hz is above median+noiseLimit*mad the piece is
% considered noisy.  the 10 s pieces are advanced in tStep and the
% computation repeats.

% Dec 2011  MA

%% initialize
if samplingRate<600
    error('MATLAB:MEGanalysis:badParam',...
        ['Sampling rate must be greater then 700 to avoid ' ...
        'confusion with line frequency artifacts'])
end
if ~exist('tStep', 'var'), tStep=[]; end
if isempty(tStep), tStep=5; end
if ~exist('lineFreq', 'var'), lineFreq=50; end
if isempty(lineFreq), lineFreq=50; end
if ~exist('noiseLimit', 'var'), noiseLimit=[]; end
if isempty(noiseLimit), noiseLimit=180; end
nChannels = size(MEG,1);
OKperiods = cell(1,nChannels);

%% Analyse
for ii = 1:nChannels
    x = MEG(ii,:);
    OKperiods{ii} = findOKperiods(x, samplingRate, noiseLimit,...
        tStep, lineFreq);
end

return
