function [wherisHB, zTime, Errors, amplitudes, meanBeat] = findHB01(mMEG, samplingRate, HBperiod, toPlot, toVerify)
% find heart beats in mean MEG channel
% [wherisHB, zTime, Errors, amplitudes, meanBeat] = findHB01(mMEG, samplingRate, HBperiod, toPlot, toVerify);
%
% mMEG         - The mean for all MEG channels
% samplingRate - of MEG data
% HBperiod     - mean interval between heart bits in s. [default 0.9];
% toPlot       - either 'plot' or 'noPlot'
% toVerify     - either 'verify' or 'noVerify'.  If 'verify' then the heart
%                beats are checked for extreme values and irregular
%                intervals.  If found the user is informed and has to
%                decide what to do. [default 'noVerify']
%
% wherisHB - list of indices for position of QRS shapes
% Errors   - struct with places of outliers
% meanBeat - shape of the mean beat
% amplitudes - list of size of QRS complexes
% zTime    - where in the cycle is time zero

% NOV 2008  MA
% UPDATES
%  Dec 2009 - tested and algorithm modified, VERIFY aded.
%  Jan 2010 - Position of zTime closer to start of cycle
%  Nov 2010 - amplitudes of initial peaks are tested and returned in Errors
%  Dec 2010 - Bug in computing HBperiod in samples - fixed
%             data is first filtered 4-50 Hz) for finding QRS position, and
%             then the mean cycle is computed on the unaltered data. 
%             test for inverted QRS. 
%             If the bigest peaks are negative the entire signal is
%             inverted.  MA

%% initialize
% internal parameters to be put as varargin in later stage
numSD=2.5;       % above that it is a QRS complex
minPeakSize=2.7; % below that it is not true qrs peak
biggest = 2;   % highest amplitude of QRS
smallest = 0.5;  % lowest amplitude of QRS

% test for missing input parameters
if ~exist('toPlot','var'), toPlot=[]; end
if isempty(toPlot), toPlot='plot'; end
if strcmpi('PLOT',toPlot)
    dbgPlot=true;  % set to TRUE if you need to see how well are the peaks detected
else
    dbgPlot=false;
end
if ~exist('toVerify','var'), toVerify=[]; end
if isempty(toVerify), toVerify='noVerify'; end
if strcmpi('VERIFY',toVerify)
    verify=true;  % set to TRUE if you need to see how well are the peaks detected
else
    verify=false;
end
Errors = struct('shortHB',[] , 'longHB',[] , 'smallHB',[] , 'bigHB',[] ...
    , 'numSmallPeaks',[]);

lastSample = size(mMEG,2);
% normalize the data
nFactor = max(mMEG);
mMEG = mMEG/nFactor;
if lastSample/samplingRate <25
    warning('MATLAB:MEGanalysis:notEnoughData','Data is too short to reliably clean the hart beat')
end
if ~exist('HBperiod','var'), HBperiod=[]; end
if isempty(HBperiod), HBperiod=1.1; end % estimate of the heart HBperiod
slowest = 1.7*HBperiod;   % lowest limit of HB period
fastest = 0.5*HBperiod;   % highest limit of HB period

beatPeriod = round(1.2*samplingRate*HBperiod);
thirdBeat = round(beatPeriod/3);
meanBeat = zeros(1,beatPeriod+1);
n=0;

%% filter the data to eliminate slow drifts and fast noise
BandPassSpecObj=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
2,4,50,100,60,1,60,samplingRate);
BandPassFilt=design(BandPassSpecObj  ,'butter');
mMEGf = myFilt(mMEG,BandPassFilt);

%% find where are the hart beats
startSample=thirdBeat;
endSample = startSample+beatPeriod;
numBeats = ceil(length(mMEG)/(thirdBeat));
pSize = nan(1,numBeats);
while endSample<lastSample-beatPeriod
    MEGpiece=mMEGf(startSample:endSample);
    peak=find(MEGpiece == max(MEGpiece),1);
    strt=startSample+peak-thirdBeat;
    meanBeat = meanBeat+mMEGf(strt:strt+beatPeriod);
    startSample = startSample+peak+thirdBeat;
    endSample = startSample+beatPeriod;
    n = n+1;
    pSize(n) = (MEGpiece(peak) - mean(MEGpiece))/std(MEGpiece);
end
Errors.numSmallPeaks = sum(pSize<=minPeakSize);
meanBeat = meanBeat/n;
maxm=find(meanBeat == max(meanBeat));
range = round(0.05*samplingRate); % 50 ms
QRS = meanBeat(maxm-range:maxm+range);
% normalize
QRS = QRS/sqrt(sum(QRS.*QRS));

%% find the position of heart beats more precisely and recompute the meanBeat
sMEG = conv(QRS,mMEGf);
sMEG = sMEG(range:end-range+1);
% we assume that QRS is pointing up!
[amplitudes, Ipeaks] = findPeaks(sMEG,numSD,thirdBeat); % find peaks >4*sd of signal
[amplitudesN, IpeaksN] = findPeaks(-sMEG,numSD,thirdBeat); % find negative peaks
if mean(amplitudesN) > mean(amplitudes) % invert
    amplitudes = amplitudesN;
    Ipeaks = IpeaksN;
    sMEG = -sMEG;
end
% check for errors
mni = min(diff(Ipeaks))/samplingRate;
if mni<0.25
    warning ('MATLAB:MEGanalysis:impossibleValue','Minimal QRS intervasl too small!')
end
mxi = max(diff(Ipeaks))/samplingRate;
if mxi>2
    warning ('MATLAB:MEGanalysis:impossibleValue','Maximal QRS intervasl too big!')
end

QRSindx = Ipeaks;
if dbgPlot
    figure
    plot(sMEG*nFactor);
    hold on
    yr = max(sMEG*nFactor);
    for ii=1:length(QRSindx)
        x=QRSindx(ii);
        y=[0.9*yr,yr];
        line([x,x],y,'color','k');
    end
    set (gcf,'Position',[6 373 1020 293]);
    set (gca,'Position',[0.04 0.174 0.95 0.72]);
    axis tight
    yl = get(gca,'YLim');
    set(gca,'YLim',[yl(1),1.1*yl(2)])
end

%% recompute the mean
meanBeat = zeros(1,ceil(samplingRate*mxi)+1);
numSamples = meanBeat;
lengthMean = length(meanBeat);
zTime = round(0.15*lengthMean); 
for ii=1:length(QRSindx)
    iStart = QRSindx(ii)-zTime;
    iEnds = QRSindx(ii) +lengthMean-zTime-1;
    numP = iEnds-iStart+1;
    if(iStart>0)&&(iEnds<=lastSample)
       meanBeat(1:numP)= meanBeat(1:numP) +mMEG(iStart:iEnds);
       numSamples(1:numP)= numSamples(1:numP)+1;
    end
end
meanBeat = nFactor*meanBeat./numSamples;
wherisHB = QRSindx;

%% verify that heart beat is regular and homogenous
if verify
    % check for timing
    dHB = diff(wherisHB);
    % HBperiod = mean(dHB);
    if max(dHB)>slowest*samplingRate
        whereSlow = wherisHB(dHB>=slowest*samplingRate)/samplingRate;
        warning('MATLAB:MEGanalysis:unreasonableValue',...
            ['Too slow Heart beat periods in ' num2str(whereSlow')]);
        Errors.longHB = whereSlow';
    end
    if min(dHB)<fastest*samplingRate
        whereFast = wherisHB(dHB<=fastest*samplingRate)/samplingRate;
        warning('MATLAB:MEGanalysis:unreasonableValue',...
            ['Too fast Heart beat periods in ' num2str(whereFast')]);
        Errors.shortHB = whereFast';
    end
    % check for amplitudes
    peaks = sMEG(wherisHB);
    meanPeak = mean(peaks);
    whereSmall = wherisHB(peaks<smallest*meanPeak)/samplingRate;
    if ~isempty(whereSmall)
        warning('MATLAB:MEGanalysis:unreasonableValue',...
            ['Too small Heart beat amplitude in ' num2str(whereSmall')]);
        Errors.smallHB = whereSmall';
    end
    whereBig = wherisHB(peaks>biggest*meanPeak)/samplingRate;
    if ~isempty(whereBig)
        warning('MATLAB:MEGanalysis:unreasonableValue',...
            ['Too large Heart beat amplitude in ' num2str(whereBig')]);
        Errors.bigHB = whereBig';
    end
    if dbgPlot  % mark the outliers
        whereMark = sort([Errors.longHB Errors.shortHB Errors.smallHB Errors.bigHB]);
        whereMark = round(samplingRate*whereMark);
        for ii=1:length(whereMark)
            x=whereMark(ii);
            y=[0.9*yr,yr];
            line([x,x],y,'color','r','LineWidth', 5);
        end
    end
        
end

%% improve the zTime by finding the center of gravity of the peak.
lHB = find(meanBeat>=meanBeat(zTime)/2,1); % half height to the left
rHB = find(meanBeat>=meanBeat(zTime)/2,1,'last'); % half height to the right
% tests show that the peak coordination is with -1
zTime = round(sum(meanBeat(lHB:rHB).*(lHB:rHB))/sum(meanBeat(lHB:rHB)))-1;

return
