function [data,HBtimes,templateHB,Period,MCG,Rtopo]=correctHB(data,sRate,figOptions,ECG,cfg)
% [data,HBtimes,templateHB,Period,MCG,Rtopo]=correctHB(data,sRate,figOptions,ECG,cfg)
% HB components
%      /\      _____          /\      _____
% __  /  \  __/     \______  /  \  __/     \__...
%   \/    \/               \/    \/
%   Q  R  S     T
%
%% Input
% - data is a matrix with rows for channels, raw data, not filtered. it can
% also be a filename.mat, directing to data matrix file, or a 4D filename such
% as 'c1,rfhp0.1Hz'.
% - sRate is sampling rate
% - figOptions=false for fuw figures (you always get them) and true for
% many figures. If you want fieldtrip topoplot of HB (first and second sweeps)
% you have to have:
% figOptions.label=N by 1 cell array with channel names of data
% figOptions.layout='4D248.lay' for 4D users.
% I recommend 'neuromag306mag.lay' for neuromag users even if data includes also grads.
% - ECG can be ECG (useful for ctf users) or a mean of subset of MEG channels where the HB is
% visible. Neuromag users can put there the mean of the magnetometers to
% clean both magnetometers and gradiometers included in data.
%
% cfg
%
% You can change variables by placing them in cfg.
% Lots of thresholds can be set and filters factors.
% a bit about the FILTERS.
% PEAK detection is performed on meanMEGf (filtered mean MEG channel). A highpass is recommended to
% ignore the T wave which is sometimes larger than the R and not very indicative of HB timing.
% cfg.peakFiltFreq sets the bandpass freqs for this one.
% The AMPLITUDE of R is estimated, but then you may want a highpass filter
% only, not to cut down R peak. cfg.ampFiltFreq can set up this filter as bp or hp.
% When using xcorr or Abeles method a TEMPLATE HB is fitted to the ECG like
% trace (meanMEG). Another filter can be used there, in order to supress T
% for better timing. When the R is small the T can help with the template
% match, so leave some low frequencies in. Use tempFiltFreq for this one.
%
%% cfg options
%  - cfg.dataFiltFreq is a filtering option  for the data (and meanMEG) to get rid of low
% frequencies (DC recordings). Not performed by default. This is the
% only filter that actually changes the data. The rest are for better
% processing of HB.
%  - cfg.chanSnrThr (default 0) is the threshold (in z scores) that tell which channels are cleaned and which
% remain as are. Use 0 to clean all. SNR here is how much heartbeat peak
% there is in each channel. You may want not to clean channels that are
% not dirty.
%  - cfg.rThr (0.5) is the threshold for correlation between topographies of averaged R peak
% and a particular R peak, it determains which instances of HB will be
% taken to calculate the HB temporal template.
%  - cfg.minPeriod (0.45) low limit for HB period. May be less for babies.
%  - cfg.maxPeriod (2) upper limit for HB period. Can be more for Mark and
%  other extreme cases, but he is now in Costa Rica
%  - cfg.chanZthr (20) z-score for declaring bad channels. A channel is bad
%  if it surpasses this value for 3 seconds in the beginning.
%  - cfg.badChan ([]) is the indices of bad channels. if you specify this, badChan
%  will not be self determined. cfg.chanZthr will be ignored. 4D users -
%  for bad A204 use cfg.badChan=204;
%  - cfg.jPad (1) how much to zero pad before and after jump
%  - cfg.jZthr (15) is a z-score over which the meanMEG is considered to
%  have a jump artifact
%  - cfg.peakFiltFreq ([7 90]) is the band-pass filter used for the meanMEG data, before peak
% detection.
%  - cfg.tempFiltFreq (same default as peakFiltFreq) is the filter used for
%  meanMEG template and meanMEG before template match takes place. When
%  T is large and R is small you may want to lower the highpass freq. It
%  can save the day but beware, tricky business.
%  - cfg.matchMethod can be 'xcorr' (default) or 'Abeles', it is how you find
% the match between a template HB and meanMEG / ECG recording. You can also
% use 'topo' and 'meanMEG' in order to define HB peaks on the topography
% trace or the mean(MEG) channel.
%  - cfg.ampFiltFreq (2) is a high-pass or band-pass filter used to test R
% peak amplitude. Better use highpass only, although bp may improve linear
% fit between (unfiltered) template and filtered data.
%  - cfg.ampLinThr (0.25) is linear regression r threshold. If there is no good
%  fit between template QRS and an instance of a heart beat, amplitude will not
%  be assesed by r, the average HB amp will be given. use cfg.ampLinThr =
%  0; for not adjusting amplitude HB by HB, remove mean template as is for
%  every HB (good when HB are rather small).
%  - cfg.ampMethod can be '5cat' (five HB size categories), '1size' to make one size template for all HB, or 'HBbyHB'
%  (correct each HB by the amplitude estimate, less recommended).
%  - cfg.afterHBs (0.7 of the period) is how long the template should continue after the
% peak (seconds)
%  - cfg.beforeHBs (0.3 of the period) is when the template should start before the
% peak (seconds)
%  - cfg.repressTime (20ms) is how much of the template to repress to zero on
%  the edges.
%
%% output
% NOTE !!! if you give no op argument MEG data will be saved
% 4D and neuromag to native format, ctf to FieldTrip structure
% (continuous).
%  - data is a matrix of cleaned data
%  - HBtimes is the times when HB peaks were detected
%  - templateHB is the average HB
%  - Period is the median interval between two heartbeas.
%  - MCG is the meanMEG channel
%  - Rtopo is the topography at the peak of R, similar to component
%  analysis topography (just cleaner)
%
%% Save File
% when no output arguments are given a 4D or a fif file will be saved!
%% Examples
% 4D users can run the function from the folder with the data ('c,*') file, with no
% input arguments:
% cleanData=correctHB;
% or like this cleanData=correctHB([],[],1); to get the figures.
% If you don't specify figure options you still get one before / after figure.
% You can write a 4D file like this: rewrite_pdf(cleanData);
%
% To clean other data use:
% cleanData=correctHB(data,samplingRate,1);
% or the more RAM friendly form
% cleanData=correctHB('data.mat',samplingRate,1);
% added by Dr. Yuval Harpaz to Prof. Abeles' work. 4D users can use Abeles'
% function createCleanFile to clean HB, but it is more touchy.


%% default variables and parameters
if ~exist('figOptions','var')
    figOptions=[];
end
if isempty(figOptions) || figOptions==0 || figOptions==false
    figs=false;
else
    figs=true;
end
if ~exist('ECG','var')
    ECG=[];
end
if ~exist('cfg','var')
    cfg=struct;
end
% uses default unless specified in the cfg
chanSnrThr=default('chanSnrThr',0,cfg);
rThr=default('rThr',0.5,cfg);
peakZthr=default('peakZthr',1.5,cfg);
minPeriod=default('minPeriod',0.45,cfg);
maxPeriod=default('maxPeriod',2,cfg);
chanZthr=default('chanZthr',8,cfg);
jPad=default('jPad',1,cfg);
jZthr=default('jZthr',15,cfg);
peakFiltFreq=default('peakFiltFreq',[7 90],cfg);
peakFiltAlpha=default('peakFiltAlpha',false,cfg);
ampFiltFreq=default('ampFiltFreq',2,cfg);% [7 90] for band pass
tempFiltFreq=default('tempFiltFreq',peakFiltFreq,cfg);
dataFiltFreq=default('dataFiltFreq',[],cfg);
matchMethod=default('matchMethod','xcorr',cfg);
beforeHBs=default('beforeHBs',[],cfg); % how long the right side of template HB should be. when empty it gets 0.3*period
afterHBs=default('afterHBs',[],cfg); % how long the right side of template HB should be. when empty it gets 0.7*period
ampLinThr=default('ampLinThr',0.25,cfg);  % threshold for low amplitude HB, use average amplitude when below this ratio
meanMEGhpFilt= default('meanMEGhpFilt',3,cfg); % highpass filter for meanMEG before everything
badChan= default('badChan',[],cfg);
ampMethod= default('ampMethod','5cat',cfg);
barilan=false; % is it Bar-Ilan University data, if so write file in the end
dataType='unknown';
%% checking defaults for 4D data
% to use with data=[] and sRate=[];
if ~exist('data','var')
    data=[];
    sRate=[];
end
if isempty(data)
    try
        rawFileName=ls('*c,rf*');
        try
            rawFileName=ls('xc,lf_c,*');
        catch %#ok<CTCH>
            try
                rawFileName=ls('lf_c,*');
            catch
                rawFileName=ls('c,rf*');
            end
        end
        data=['./',rawFileName(1:end-1)];
        dataType='4D';
    catch
        try
            rawFileName=ls('*raw.fif');
            if length(regexp(rawFileName,'raw.fif'))==1
                data=['./',rawFileName(1:end-1)];
                dataType='fif';
            else
                fprintf(2,'more than one raw.fif files? specify the name please thank you very much')
                return
            end
        end
        PWD=pwd;
        if strcmp(PWD(end-2:end),'.ds')
            data=PWD;
        else
            ctfDir=dir('*.ds');
            if length(ctfDir)~=1
                error('cannot figure out which data type it is, give me the data matrix or filename')
            else
                cd(ctfDir.name)
                data=pwd;
                cd ../
            end
        end
    end
end
if ischar(data)
    if strcmp(data(end-2:end),'fif')
        dataType='fif';
    elseif strcmp(data(end-2:end),'.ds')
        dataType='ctf';
    else
        if ~strcmp(dataType,'4D') && ~strcmp(data(end-3:end),'.mat')
            if ~isempty(findstr(',rf',data))
                dataType='4D';
            end
        end
    end
    % can be .mat or raw MEG data filename
    if strcmp(data(end-3:end),'.mat') % read matrix from file 'data.mat'
        PWD=pwd;
        display(['loading ',PWD,'/',data,]);
        data=load(['./',data]);
        dataField=fieldnames(data);
        
        eval(['data=data.',dataField{1,1},';']);
    elseif strcmp(dataType,'4D') % read 4D data from file name specified in 'data'
        var4DfileName=data;
        ipFileName=var4DfileName;
        doX=false;
        var4Dp=pdf4D(var4DfileName);
        sRate=double(get(var4Dp,'dr'));
        var4Dhdr = get(var4Dp, 'header');
        var4DnSamp=var4Dhdr.epoch_data{1,1}.pts_in_epoch;
        var4Dchi = channel_index(var4Dp, 'meg', 'name');
        cnf = get(var4Dp, 'config');
        if strcmp(cnf.config_data.site_name,'Bar Ilan')
            barilan=true;
            var4DchiX = channel_index(var4Dp, 'EXTERNAL', 'name');
            lastMEG=length(var4Dchi);
            if length(var4DchiX)>5
                doX=true;
                var4Dchi=[var4Dchi,var4DchiX(4:6)];
                labels=channel_name(var4Dp, var4Dchi);
                %data(end+1:end+3,:) = read_data_block(var4Dp,[1 var4DnSamp],var4Dchi);
            end
        end
        display(['reading ',var4DfileName]);
        data = read_data_block(var4Dp,[1 var4DnSamp],var4Dchi);
        if ischar(ECG)
            var4Dchi = channel_index(var4Dp, ECG, 'name');
            ECG = read_data_block(var4Dp,[1 var4DnSamp],var4Dchi);
        end
        %data=double(data);
        if figs
            var4Dlabel=channel_label(var4Dp,var4Dchi)';
            if figOptions==1;
                clear figOptions
            end
            figOptions.label=var4Dlabel;
            figOptions.layout='4D248.lay';
        end
        clear var4D*
    elseif strcmp(dataType,'fif')
        if isempty(which('fiff_read_raw_segment'))
            error('trying to read neuromag data requires fiff_read_raw_segment.m in path');
        else
            infile=data;
            if strcmp(infile(1:2),'./')
                outfile=['./hb,',infile(3:end)];
            else
                outfile=['hb,',infile];
            end
            raw=fiff_setup_read_raw(data);
            [data,times] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp);
            sRate=1/diff(times(1:2));
            clear times
            magCount=0;
            megCount=0;
            for chani=1:length(raw.info.ch_names)
                if strcmp(raw.info.ch_names{chani}(1:3),'MEG')
                    megCount=megCount+1;
                    megi(megCount)=chani;
                    if strcmp(raw.info.ch_names{chani}(end),'1')
                        magCount=magCount+1;
                        magi(magCount)=chani;
                    end
                end
            end
            ECG=mean(data(magi,:));
            data=data(megi,:);
        end
    elseif strcmp(dataType,'ctf')
        if isempty(which('ctf_read'))
            error('to read ctf data you should have ctf_read.m in path')
        else
            ctf = ctf_read(data);
            megi=ctf.sensor.index.meg;
            data=ctf.data{1}(:,megi)';
            if length(ctf.data)>1
                for segi=2:length(ctf.data)
                    data=[data,ctf.data{segi}(:,megi)'];
                end
            end
            if isempty(find(data(:,end)))
                realEnd=find(fliplr(sum(data)),1);
                data=data(:,1:size(data,2)-realEnd+1);
            end
            sRate=ctf.setup.sample_rate;
            ecgi=[];
            % try reading ECG data
            if isempty(ECG)
                [~,ecgi]=ismember('ECG',ctf.sensor.label);
                if ecgi==0
                    ecgi=ctf.sensor.index.eeg;
                end
            elseif ischar(ECG)
                [~,ecgi]=ismember(ECG,ctf.sensor.label);
            elseif length(ECG)==1
                ecgi=ECG;
            end
            
            if ~isempty(ecgi) % ECG channel not supplied and some EEG channels found
                ecg=ctf.data{1}(:,ecgi)';
                if length(ctf.data)>1
                    for segi=2:length(ctf.data)
                        ecg=[ecg,ctf.data{segi}(:,ecgi)'];
                    end
                end
                if exist('realEnd','var')
                    ecg=ecg(:,1:size(ecg,2)-realEnd+1);
                end
                if size(ecg,1)>1
                    %kur=kurtosis(ecg');
                    kur=kurLoop(ecg,sRate);
                    [~,maxKur]=max(kur);
                    ecg=ecg(maxKur,:);
                end
                ECG=ecg;
                clear ecg
            end
            ctf=rmfield(ctf,'data');
            %clear ctf
            cd(PWD)
        end
    end
end
if ~exist('lastMEG','var')
    lastMEG=size(data,1);
end

if ~isempty(ECG)
    fECG=abs(fftBasic(ECG,sRate));
    blECG=mean(fECG(5:25));
    if fECG(60)/blECG>2 || fECG(50)/blECG>2
        ECG=correctLF(ECG,sRate);
        title('cleaned ECG channel')
    end
    meanMEG=ECG;
    %meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
else
    meanMEG=double(mean(data(1:lastMEG,:)));
end
meanMEG=meanMEG-median(meanMEG(1:round(sRate)));
if ~isempty(meanMEGhpFilt)
    mmhpObj  = fdesign.highpass('N,F3dB', 10, meanMEGhpFilt, sRate);
    mmhpFilt=design(mmhpObj ,'butter');
    meanMEG=myFilt(meanMEG,mmhpFilt);
end
repressTime=default('repressTime',round(sRate/50),cfg);
%% pad with zeros (or with some baseline constant) for template slide
sampBefore=round(sRate*maxPeriod);
time=-sampBefore/sRate:1/sRate:(size(data,2)+sampBefore)/sRate;
time=time(2:end);
BL0=median(meanMEG(1:sampBefore));
BL1=median(meanMEG(end-sampBefore+1:end));
meanMEG=[ones(1,sampBefore)*BL0,meanMEG,ones(1,sampBefore)*BL1];
%meanMEG(end-sampBefore+1:end)=median(meanMEG(end-2*sampBefore:end-sampBefore));
BL0=median(data(:,1:sampBefore),2);
BL1=median(data(:,end-sampBefore+1:end),2);
data=[repmat(BL0,1,sampBefore),data,repmat(BL1,1,sampBefore)];
%data(:,end-sampBefore+1:end)=repmat(median(data(:,end-2*sampBefore:end-sampBefore),2),1,sampBefore);
%% Filter data if requested
if ~isempty(dataFiltFreq)
    display('filtering the data')
    if length(dataFiltFreq)==2
        Fp1=dataFiltFreq(1);
        Fst1=max([0.01,Fp1-1]);
        ObjData=fdesign.bandpass(...
            'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
            Fst1,Fp1,dataFiltFreq(2),dataFiltFreq(2)+10,60,1,60,sRate);
    elseif length(dataFiltFreq)==1
        if dataFiltFreq<15
            ObjData=fdesign.highpass('Fst,Fp,Ast,Ap',max([dataFiltFreq-1,0.01]),dataFiltFreq,60,1,sRate);%
        elseif dataFiltFreq>=15
            ObjData=fdesign.lowpass('Fp,Fst,Ap,Ast',dataFiltFreq,dataFiltFreq+10,1,60,sRate);
        end
    end
    FiltData=design(ObjData ,'butter');
    data = myFilt(data,FiltData);
    meanMEG=myFilt(meanMEG,FiltData);
end

for chani=1:size(data,1)
    data(chani,:)=data(chani,:)-median(data(chani,:));
end
meanMEG=meanMEG-median(meanMEG);
%% filter mean MEG (or ECG)
% filtering to pass from 5-7Hz to 90-110Hz
if peakFiltAlpha
    %bsObj = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2',1,9,11,100,1,60,1,sRate);
    %bsObj = fdesign.bandstop('N,Fp1,Fst1,Fst2,Fp2',9,1,9,11,100,sRate);
    %bsFilt = design(bsObj,'butter');
    bsObj  = fdesign.bandstop('N,F3dB1,F3dB2', 10, 8, 12, sRate);
    bsFilt = design(bsObj, 'butter');
    meanMEGf = myFilt(meanMEG,bsFilt);
end
if ~(peakFiltFreq(1)>1)
    Fst1=0.001;
else
    Fst1=peakFiltFreq(1)-1;
end
% BandPassSpecObj=fdesign.bandpass(...
%     'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
%     Fst1,peakFiltFreq(1),peakFiltFreq(2),peakFiltFreq(2)+10,60,1,60,sRate);
BandPassSpecObj  = fdesign.bandpass('N,F3dB1,F3dB2', 10, peakFiltFreq(1), peakFiltFreq(2), sRate);
BandPassFilt=design(BandPassSpecObj ,'butter');
if peakFiltAlpha
    meanMEGf = myFilt(meanMEGf,BandPassFilt);
else
    meanMEGf = myFilt(meanMEG,BandPassFilt);
end
% baseline correction again, just in case
meanMEGf=meanMEGf-median(meanMEGf);


%% look for a noisy segments and noisy channels
% find bad channels, has to be noisy for 3 of the first 3 seconds
if isempty(badChan)
    stdMEG=std(data(1:lastMEG,1+sampBefore:round(sRate)+sampBefore)');
    badc1=stdMEG>chanZthr*median(stdMEG);
    stdMEG=std(data(1:lastMEG,round(sRate)+sampBefore:round(2*sRate)+sampBefore)');
    badc2=stdMEG>chanZthr*median(stdMEG);
    stdMEG=std(data(1:lastMEG,round(2*sRate)+sampBefore:round(3*sRate)+sampBefore)');
    badc3=stdMEG>chanZthr*median(stdMEG);
    badc=find((badc1+badc2+badc3)==3);
else
    badc=badChan;
end
if ~isempty(badc)
    data(badc,:)=0;
end

% find jump or other huge artifact
zMEG=(meanMEGf-mean(meanMEGf))./std(meanMEGf);
jbeg=find(abs((zMEG))>jZthr,1);%
j=find(abs((zMEG))>jZthr);
bads=[]; % bad samples
if ~isempty(jbeg)
    for jumpi=1:length(j)
        bads=[bads,j(jumpi)-round(sRate*jPad):j(jumpi)+sRate*jPad]; %#ok<AGROW>
    end
    bads=unique(bads);
    bads=bads(bads>0);
    badData=data(:,bads);
    if unique(diff(bads))==1
        % here is an attempt to avoid hp filter, setting bad data not to
        % zero but to a line connecting the good parts. works only for one
        % bad segment
        for chani=1:size(data,1)
            if ~ismember(chani,badc)
                BL0=mean(data(chani,(bads(1)-round(sRate*0.1)):(bads(1)-1)));
                BL1=mean(data(chani,(bads(end)+1):(bads(end)+round(sRate*0.1))));
                linConnect=BL0:((BL1-BL0)/length(bads)):BL1;
                data(chani,bads)=linConnect(2:end);
            end
        end
    else
        data(chani,bads)=0;
    end
    if length(data)<2^19
        data=data-repmat(median(data,2),1,size(data,2));
    else
        for chani=1:size(data,1)
            data(chani,:)=data(chani,:)-median(data(chani,:),2);
        end
    end
    diary('HBlog.txt')
    warning(['jump? what''','s the noise at ',num2str(max([0,time(bads(1))])),'s? zeroed noise from ',num2str(max([0,time(bads(1))])),' to ',num2str(time(bads(end)))]);
    diary off
    % end
else
    % baseline correction by removing the median of each channel
    if length(data)<2^19
        data=data-repmat(median(data,2),1,size(data,2));
    else
        for chani=1:size(data,1)
            data(chani,:)=data(chani,:)-median(data(chani,:),2);
        end
    end
end
if isempty(ECG)
    meanMEG=double(mean(data(1:lastMEG,:)));
    
    meanMEG(end-sampBefore+1:end)=median(meanMEG(end-2*sampBefore:end-sampBefore));
    if ~isempty(meanMEGhpFilt)
        meanMEG=meanMEG-mean(meanMEG(1:sampBefore));
        meanMEG=myFilt(meanMEG,mmhpFilt);
    end
    meanMEGf = myFilt(meanMEG,BandPassFilt);
    meanMEGf=meanMEGf-median(meanMEGf);
end
%% peak detection on MCG (or ECG) signal
disp('looking for HB peaks')
[peaks, Ipeaks]=findPeaks(meanMEGf,peakZthr,round(sRate*minPeriod)); % 450ms interval minimum
% test if, by chance, the HB field is mainly negative
posHB=true;
[peaksNeg, IpeaksNeg]=findPeaks(-meanMEGf,peakZthr,round(sRate*minPeriod));
if median(peaksNeg)/median(peaks)>1.1
    diary('HBlog.txt')
    warning('NEGATIVE HEARTBEATS?');
    diary off
    period1=median(diff(IpeaksNeg))./sRate;
    if period1<2
        [peaks, Ipeaks]=findPeaks(-meanMEGf,peakZthr,round(sRate*period1*0.6));
        peaks=-peaks;
    else
        Ipeaks=IpeaksNeg;
        peaks=-peaksNeg;
    end
    posHB=false;
    %meanMEGf=-meanMEGf;
else
    period1=median(diff(Ipeaks))./sRate;
    if period1<2 %#ok<*BDSCI> %try to improve peak detection if peak intervals are reasonable
        [peaks, Ipeaks]=findPeaks(meanMEGf,peakZthr,round(sRate*period1*0.6)); % 450ms interval
    end
end

if figs
    figure;
    plot(time,meanMEGf)
    hold on
    plot(time(Ipeaks), peaks,'ro')
    title('peak detection on mean MEG (or ECG) trace, OK if many of them are not HB')
end
%% get topography
% if figs
%     if isfield(figOptions,'layout') && isfield(figOptions,'label')
%         topo={};
%         topo.avg=median(data(1:lastMEG,Ipeaks),2);
%         topo.time=0;
%         topo.label=figOptions.label;
%         topo.dimord='chan_time';
%         cfgp=[];
%         cfgp.layout=figOptions.layout;
%         if ~isempty(badc)
%             cfgp.channel=setdiff(1:length(topo.label),badc);
%         end
%         if strcmp(cfgp.layout,'neuromag306mag.lay')
%             [~,magi]=ismember({'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641'},topo.label);
%             %topo.avg=topo.avg(chi);
%             %topo.label=topo.label(chi);
%             cfgp.xlim=[1,1];
%             cfgp.zlim=[-max(abs(topo.avg(magi))) max(abs(topo.avg(magi)))];
%             figure;
%             ft_topoplotER(cfgp,topo);
%             title('MAGNETOMETERS, TOPOGRAPHY OF R')
%             % cfg.layout='neuromag306planar.lay';
%             % grd=topo.avg;
%             % grd(chi)=0;
%             % cfg.zlim=[-max(abs(grd)) max(abs(grd))];
%             % figure;
%             % ft_topoplotER(cfg,topo);
%             % title('GRADIOMETERS')
%         else
%             %cfg.channel={'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641'};
%             %cfg.interpolation='linear';
%             cfgp.xlim=[1,1];
%             cfgp.zlim=[-max(abs(topo.avg)) max(abs(topo.avg))];
%             figure;
%             ft_topoplotER(cfgp,topo);
%             title('TOPOGRAPHY OF R')
%         end
%     else
%         warning('no topoplot without labels and layout fields! see figOptions options')
%     end
% end
topoTrace=median(data(1:lastMEG,Ipeaks),2)'*data(1:lastMEG,:);
topoTrace(end-sampBefore+1:end)=median(topoTrace(end-2*sampBefore:end-sampBefore));
topoTrace=topoTrace-mean(topoTrace(1:sampBefore));
topoTrace=myFilt(topoTrace,BandPassFilt);
topoTrace=topoTrace-median(topoTrace);
% if ~posHB
%     meanMEGN=meanMEGf./max(-meanMEGf(1:round(sRate*10)));
%     topoTraceN=-topoTrace./max(topoTrace(1:round(sRate*10)));
% else
%     topoTraceN=topoTrace./max(topoTrace(1:round(sRate*10)));
%     meanMEGN=meanMEGf./max(meanMEGf(1:round(sRate*10)));
% end
meanMEGN=standard(meanMEGf);
topoTraceN=standard(topoTrace);
if ~posHB
    topoTraceN=-topoTraceN;
end


%% check if topo of every peak is correlated to average topo
r=corr(data(1:lastMEG,Ipeaks),median(data(1:lastMEG,Ipeaks),2));
if figs
    figure;
    plot(time,topoTraceN)
    hold on
    plot(time,meanMEGN,'r')
    plot(time(Ipeaks(r>rThr)),meanMEGN(Ipeaks(r>rThr)),'g.');
    legend('topoTrace','meanMEG',['r data-topo > ',num2str(rThr)])
end
%% average good HB and make a template
IpeaksR=Ipeaks(r>0.5);
if length(IpeaksR)<(length(Ipeaks/2))
    IpeaksR=Ipeaks(r>prctile(r,25));
    diary('HBlog.txt')
    disp(['poor correlation between median topography and each HB, adjusting rThr to ',num2str(prctile(r,25))]);
    diary off
end
IpeaksR=IpeaksR(IpeaksR>sampBefore);
IpeaksR=IpeaksR(IpeaksR<(size(data,2)-sampBefore));
period2=diff(IpeaksR)/sRate; % estimate period
period2=median(period2(period2<maxPeriod)); % less than 2s

%LowPassSpecObj=fdesign.lowpass('Fp,Fst,Ap,Ast',45,55,1,60,sRate);
if tempFiltFreq==peakFiltFreq
    meanMEGxcrF=meanMEGf;
else
    if length(tempFiltFreq)==2
        Fp1=tempFiltFreq(1);
        Fst1=max([0.01,Fp1-1]);
        ObjXcr=fdesign.bandpass(...
            'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
            Fst1,Fp1,tempFiltFreq(2),tempFiltFreq(2)+10,60,1,60,sRate);
        FiltXcr=design(ObjXcr ,'butter');
        meanMEGxcrF = myFilt(meanMEG,FiltXcr);
    elseif length(tempFiltFreq)==1
        if tempFiltFreq<15
            ObjXcr=fdesign.highpass('Fst,Fp,Ast,Ap',max([tempFiltFreq-1,0.01]),tempFiltFreq,60,1,sRate);%
        elseif tempFiltFreq>=15
            ObjXcr=fdesign.lowpass('Fp,Fst,Ap,Ast',tempFiltFreq,tempFiltFreq+10,1,60,sRate);
        end
        FiltXcr=design(ObjXcr ,'butter');
        meanMEGxcrF = myFilt(meanMEG,FiltXcr);
    end
end


% meanMEGxcrF = myFilt(meanMEG,BandPassFiltXcr);
meanMEGxcrF=meanMEGxcrF-median(meanMEGxcrF);
if posHB
    [temp1e,period3]=makeTempHB(meanMEGxcrF,sRate,IpeaksR,period2,sampBefore,figs,maxPeriod);
else
    [temp1e,period3]=makeTempHB(-meanMEGxcrF,sRate,IpeaksR,period2,sampBefore,figs,maxPeriod);
end
%% find xcorr between template and meanMEG
maxi=round(0.3*length(temp1e))+1; %maxi is where the R peak in the template
[~,mi]=max(temp1e(maxi-round(sRate/100):maxi+round(sRate/100)));
maxi=maxi-round(sRate/100)+mi-1;
if posHB
    meanMEGpos=meanMEGxcrF;
else
    meanMEGpos=-meanMEGxcrF;
end
switch matchMethod
    case 'xcorr'
        matchTrace=XCORR(meanMEGpos,temp1e);
    case 'Abeles'
        [snr,~]=match_temp(meanMEGxcrF,temp1e,maxi);
        matchTrace=snr;
    case 'topo'
        matchTrace=topoTrace;
    case 'meanMEG'
        matchTrace=meanMEGf;
end


if figs
    figure;
    plot(time,topoTraceN)
    hold on
    plot(time,meanMEGN,'r')
    if posHB
        plot(time,standard(matchTrace),'g');
    else
        plot(time,-standard(matchTrace),'g');
    end
    legend('topoTrace','meanMEG','temp xcorr')
end
%% second sweep
%% find peaks on xcorr trace
[peaks2, Ipeaks2]=findPeaks(matchTrace,peakZthr,round(sRate*period3*0.65)); % no peaks closer than 60% of period
if figs
    figure;
    if posHB
        plot(time,matchTrace)
        hold on
        plot(time(Ipeaks2), peaks2,'ro')
    else
        plot(time,-matchTrace)
        hold on
        plot(time(Ipeaks2), -peaks2,'ro')
    end
    title('2nd sweep peak detection, based on template matching')
end
% checking results
periodS=diff(Ipeaks2);% the period in samples for each HB
farI=find(periodS/median(periodS)>1.5);
if median(periodS)/sRate<0.5
    diary('HBlog.txt')
    warning('interval between HB is less than 0.5s, look at the figures!')
    diary off
elseif median(periodS)/sRate>1.5
    diary('HBlog.txt')
    warning('interval between HB is more than 1.5s, look at the figures!')
    diary off
else
    if ~isempty(farI)
        farT=Ipeaks2(farI)/sRate; % this is far time, not what you think!
        diary('HBlog.txt')
        disp('sparse heartbeats, you may want to look for missing HB after:');
        disp(num2str(farT));
        diary off
    end
    nearI=find(periodS/median(periodS)<0.5);
    if ~isempty(nearI)
        nearT=Ipeaks2(nearI)/sRate;
        diary('HBlog.txt')
        disp('close heartbeats, you may want to look for false HB detection after:');
        disp(num2str(nearT));
        diary off
    end
end
% ignore edges
Ipeaks2in=Ipeaks2(Ipeaks2>sampBefore);
Ipeaks2in=Ipeaks2in(Ipeaks2in<(size(data,2)-sampBefore));
% set template edges as ratio of the period
if isempty(afterHBs)
    afterHBs=0.7;
else
    afterHBs=afterHBs/period3;
end
if isempty(beforeHBs)
    beforeHBs=0.3;
else
    beforeHBs=beforeHBs/period3;
end
% make mcg trace for meanMEG
if posHB
    [templateHB,Period]=makeTempHB(meanMEG,sRate,Ipeaks2in,period3,sampBefore,figs,maxPeriod,beforeHBs,afterHBs,repressTime);
else
    [templateHB,Period]=makeTempHB(-meanMEG,sRate,Ipeaks2in,period3,sampBefore,figs,maxPeriod,beforeHBs,afterHBs,repressTime);
end

%maxi=round(0.3*length(templateHB))+1;
maxi=round(beforeHBs*Period*sRate)+1;
[~,mi]=max(templateHB(maxi-round(sRate/100):maxi+round(sRate/100)));
maxi=maxi-round(sRate/100)+mi-1;
%% test R amplitude

% %make temp per chan
% sBef=maxi-1;
% sAft=length(templateHB)-maxi;
% HBtemp=HBbyChan(data,Ipeaks2in,sBef,sAft,repressTime);
% meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
%
if length(ampFiltFreq)==1
    HighPassSpecObj=fdesign.highpass('Fst,Fp,Ast,Ap',ampFiltFreq-1,ampFiltFreq,60,1,sRate);%
    HighPassFilt=design(HighPassSpecObj ,'butter');
    meanMEGampF = myFilt(meanMEG,HighPassFilt);
    % baseline correction again, just in case
    meanMEGampF=meanMEGampF-median(meanMEGampF);
elseif length(ampFiltFreq)==2
    if isequal(ampFiltFreq,[7 90])
        meanMEGampF=meanMEGf;
    else
        BandPassSpecObjAmp=fdesign.bandpass(...
            'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
            ampFiltFreq(1)-1,ampFiltFreq(1),ampFiltFreq(2),ampFiltFreq(2)+10,60,1,60,sRate);
        BandPassFiltAmp=design(BandPassSpecObjAmp ,'butter');
        meanMEGampF = myFilt(meanMEG,BandPassFiltAmp);
        % baseline correction again, just in case
        meanMEGampF=meanMEGampF-median(meanMEGampF);
    end
elseif isempty(ampFiltFreq)
    meanMEGampF=meanMEG;
    disp('not filtering meanMEG for testing R amplitude!')
else
    error('wrong length of vector. one number (2) means hp filter, two ([2 90]) means bp')
end

[p,Rlims]=assessAmp(templateHB,maxi,Ipeaks2in,meanMEGampF);
% % if strcmp(ampMethod,'HBbyHBtopo')
% % [templateHBtopo,PeriodTopo]=makeTempHB(meanMEGampF,sRate,Ipeaks2in,period3,sampBefore,figs,maxPeriod,beforeHBs,afterHBs,repressTime);
% % [p,Rlims]=assessAmp(templateHBtopo,maxi,Ipeaks2in,meanMEGampF);

% look for neg correlation between template peak and peaks, means trouble
% if there are too many
if ampLinThr==0
    ampMMfit=ones(size(p,1),1);
else
    if posHB
        negp=find(p(:,1)<ampLinThr);
    else
        negp=find(p(:,1)>ampLinThr);
    end
    if ~isempty(negp)
        p(negp,1:2)=0;
        diary('HBlog.txt')
        warning(['did not get good fit for amplitude test, assume average HB amplitude at ',num2str(time(Ipeaks2in(negp)))])
        diary off
    end
    ampMMfit=p(:,1)+(p(:,2)./templateHB(maxi));
    ampMMfit(ampMMfit>1.5)=1;
end


if posHB
    MCG=makeMCG(templateHB,maxi,Rlims,Ipeaks2in,ampMMfit,length(meanMEG));
else
    MCG=makeMCG(templateHB,maxi,Rlims,Ipeaks2in,-ampMMfit,length(meanMEG));
    MCG=-MCG;
end

%% remove mcg from each chan
%make temp per chan
sBef=maxi-1;
sAft=length(templateHB)-maxi;
HBtemp=HBbyChan(data,Ipeaks2in,sBef,sAft,repressTime);
% check SNR per chan
sm50=round(sRate*0.05);
s0=maxi-sm50*3;
s1=maxi-sm50;
s2=maxi+sm50;
n=std(HBtemp(:,s0:s1)'); %#ok<*UDIM>
s=std(HBtemp(:,s1:s2)');
snr=s./n;
badSNR=find(snr<=chanSnrThr);
if figs
    timeTemp=1/sRate:1/sRate:size(HBtemp,2)/sRate;
    timeTemp=timeTemp-maxi/sRate;
    figure;plot(timeTemp,HBtemp','k');hold on;
    %plot(maxi,HBtemp(find(HBsnr>1.1),maxi),'g.')
    [~,minI]=min(snr);
    plot(timeTemp,HBtemp(minI,:),'r');
    if isempty(badSNR) || chanSnrThr==0
        title('HB for all channels,red trace is worst SNR channel')
    else
        plot(timeTemp(maxi),HBtemp(badSNR,maxi),'r.')
        title(['HB for all channels. red trace is worst SNR channel, dots mark channels with SNR < ',num2str(chanSnrThr)])
    end
end
if ~isempty(badSNR)
    if length(find(snr<=chanSnrThr))>size(data,1)/2
        error('too many channels have poor SNR. was there artifact?')
    end
    if ~chanSnrThr==0
        HBtemp(badSNR,:)=0;
        diary('HBlog.txt')
        disp(['not including bad SNR channels, index no. ',num2str(badSNR)])
        diary off
    end
end
% clear some memory
meanData=mean(data(1:lastMEG,:));
MEGmean=meanMEG;
clear meanMEG* topoTra*
meanMEG=MEGmean;
clear MEGmean;
% prepare avg HB fig
HBtimes=(Ipeaks2in-sampBefore)/sRate;
[avgHB,avgTimes]=meanHB(data(1:lastMEG,sampBefore+1:end-sampBefore),sRate,HBtimes);
Rtopo=HBtemp(:,maxi);
sign=2*((Rtopo>0)-0.5);
if strcmp(dataType,'ctf')
    meanDataSigned=sign'*data;
end
% clean
display('cleaning channels from template one by one, may take half a minute')
if max(max(templateHB))/max(max(avgHB))>10000 % ECG template MEG signal
    if ampLinThr~=0
        warning('will not try to fit amplitude of each HB, I think the template is ECG and the data is MEG')
        ampLinThr=0;
    end
end
if strcmp(ampMethod,'1size')
    ampLinThr=0;
end
if ampLinThr==0
    ampMethod='HBbyHB';
    disp('Using one size templae for all heartbeats')
    %not to make 5 categories
end
switch ampMethod
    case '5cat' % 5 categories of HB
        [~,sorted]=sort(ampMMfit);
        len=floor(length(sorted)/5);
        for chani=1:size(data,1)
            MCGall=zeros(size(meanMEG));
            for cati=1:5
                if cati<5
                    HBcat=sorted(len*(cati-1)+1:len*cati);
                else
                    HBcat=sorted(len*(cati-1)+1:end);
                end
                HBtempCat=HBbyChan(data(chani,:),Ipeaks2in(HBcat),sBef,sAft,repressTime);
                %HBtempCat=HBbyChan(data(1:lastMEG,:),Ipeaks2in(HBcat(:,cati)),sBef,sAft,repressTime);
                HBtempCat(1:Rlims(1))=HBtemp(chani,1:Rlims(1));
                HBtempCat(Rlims(2):end)=HBtemp(chani,Rlims(2):end);
                HBtempCatAll(chani,1:size(HBtempCat,2),cati)=HBtempCat;
                amp1=ones(length(HBcat),1);
                
                if posHB
                    MCGall=MCGall+makeMCGbyCh(HBtempCat,maxi,Rlims,Ipeaks2in(HBcat),amp1,length(meanMEG));
                else
                    MCGall=MCGall+makeMCGbyCh(HBtempCat,maxi,Rlims,Ipeaks2in(HBcat),-amp1,length(meanMEG));
                end
                
            end
            data(chani,:)=data(chani,:)-MCGall;
        end
        if figs
            figure;
            plot(squeeze(mean(HBtempCatAll)));
            legend(num2str([1:5]'));
            title('average HB by size category')
        end
    case 'HBbyHB'
        
        for chani=1:size(HBtemp,1)
            if posHB
                MCGall=makeMCGbyCh(HBtemp(chani,:),maxi,Rlims,Ipeaks2in,ampMMfit,length(meanMEG));
            else
                MCGall=makeMCGbyCh(HBtemp(chani,:),maxi,Rlims,Ipeaks2in,-ampMMfit,length(meanMEG));
            end
            data(chani,:)=data(chani,:)-MCGall;
        end
end
%cleanData=data-MCGall;

figure;

if strcmp(dataType,'ctf')

    scale=max(abs(MCG(sampBefore+1:sampBefore+round(sRate*5))))/max(abs(meanDataSigned(sampBefore+1:sampBefore+round(sRate*5))));
    plot(time,meanMEG/scale,'k')
    hold on
    plot(time,meanDataSigned,'r')
    plot(time,sign'*data,'g');
else
    if isempty(ECG)
        plot(time,MCG,'k')
    else
        scale=max(abs(MCG(sampBefore+1:sampBefore+round(sRate*5))))/max(abs(mean(data(1:lastMEG,sampBefore+1:sampBefore+round(sRate*5)))));
        plot(time,meanMEG/scale,'k')
    end
    hold on
    plot(time,meanData,'r')
    plot(time,mean(data(1:lastMEG,:)),'g')
end
if isempty(ECG)
    legend('MCG from template', 'mean MEG','mean clean MEG')
else
    legend('rescaled ECG', 'mean MEG','mean clean MEG')
end

if figs
    if isfield(figOptions,'layout')
        topo.avg=Rtopo;
        cfgp.xlim=[1,1];
        if strcmp(cfgp.layout,'neuromag306mag.lay')
            cfgp.zlim=[-max(abs(topo.avg(magi))) max(abs(topo.avg(magi)))];
            figure;
            ft_topoplotER(cfgp,topo);
            title('MAGNETOMETERS, TOPOGRAPHY OF R, 2nd sweep')
        else
            cfgp.zlim=[-max(abs(Rtopo)) max(abs(Rtopo))];
            figure;
            ft_topoplotER(cfgp,topo);
            title ('TOPOGRAPHY OF R, 2nd sweep')
        end
    end
end
avgHBclean=meanHB(data(1:lastMEG,sampBefore+1:end-sampBefore),sRate,HBtimes);
if strcmp(dataType,'ctf')
    avgHBclean=meanHB(data(1:lastMEG,sampBefore+1:end-sampBefore),sRate,HBtimes);
    avgHBclean=sign'*avgHBclean;
    figure;plot(avgTimes,sign'*avgHB,'r')
    hold on
    plot(avgTimes,avgHBclean,'g')
    title('averaged heartbeat, before (red) and after')
else
    figure;plot(avgTimes,mean(avgHB),'r')
    hold on
    plot(avgTimes,mean(avgHBclean),'g')
    title('averaged heartbeat, before (red) and after')
end
%display(['HB period (2nd sweep) is ',num2str(period4),'s']);
if ~isempty(bads);
    data(:,bads)=badData;
end
% plot avg HB before and after

data=data(:,sampBefore+1:end-sampBefore);
%% sanity check
[~,t0]=max(diff(smooth(templateHB,20)));
if t0>round(sRate/2)
    diary('HBlog.txt')
    txt=['insane, max(diff(templateHB)) is far, ',num2str(round(1000*t0/sRate)),'ms. Is it HeartBeat at all?!?\n'];
    fprintf(2,txt);
    diary('off')
end
EK=kurtosis(xcorr(templateHB))-3; %excess kurtosis of autocorrelation should be high
if EK <3
    txt='insane, repetitive pattern, is it HeartBeat? not Alpha?!?\n';
    diary('HBlog.txt')
    fprintf(2,txt);
    diary off
end
if nargout==0
    switch dataType
        case '4D'
            display('saving hb* file')
            if doX
                rewrite_pdf(data,labels,ipFileName,'hb')
            else
                rewrite_pdf(data,[],ipFileName,'hb')
            end
        case 'fif'
            copyfile(infile,outfile);
            [outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
            dataAll = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp);
            dataAll(megi,:)=data;
            clear data
            fiff_write_raw_buffer(outfid,dataAll,cals);
            fiff_finish_writing_raw(outfid);
            data=dataAll(megi,:);
            clear dataAll
        case 'ctf'
            prepPath=which('ft_preprocessing');
            if isempty(prepPath)
                warning('if you had ft_preprocessing function in path I would have saved the data to a fieldtrip structure')
            else
                cfg=[];
                cfg.dataset=pwd;
                cfg.channel=ctf.sensor.label{1};
                cfg.continuous='yes';
                cfg.feedback='none';
                cfg.trl=[1 5 1];
                correctHBdata=ft_preprocessing(cfg);
                correctHBdata.trial={};
                correctHBdata.trial{1}=data;
                clear data
                correctHBdata.label=ctf.sensor.label(megi)';
                correctHBdata.time={};
                correctHBdata.time{1}=0:1/sRate:size(correctHBdata.trial{1},2)./sRate-1/sRate;
                display('saving correctHBdata.mat, FieldTrip structure with the cleaned stuff in it')
                save correctHBdata correctHBdata ECG -v7.3
            end
    end
end

%% internal functions
function [tempe,period4]=makeTempHB(trace,sRate,peakIndex,period,sampBefore,figs,maxPeriod,beforeHBs,afterHBs,repressTime)
if ~exist('afterHBs','var')
    afterHBs=0.7; % after the T wave, before next qrs, 0.7 of period
end
if ~exist('beforeHBs','var')
    beforeHBs=0.3;
end
if ~exist('repressTime','var')
    repressTime=20; % how much time to supress
end
HB=zeros(size(trace,1),sampBefore*2+1);
for epochi=1:length(peakIndex)
    HB=HB+trace(:,peakIndex(epochi)-sampBefore:peakIndex(epochi)+sampBefore);
end
HB=HB/epochi;
period4=[]; %#ok<NASGU>
HBxc=xcorr(HB);
[~,ipxc]=findPeaks(HBxc,1.5,sRate*period*0.6); % index peak xcorr
if length(ipxc)>1
    nextPeak=ceil(length(ipxc)/2+0.5);
    period4=(ipxc(nextPeak)-ipxc(nextPeak-1))/sRate; %#ok<NASGU>
    % else
    xcorrByseg=false;
    if length(trace)>=2^20 % test if version is later than 1011b
        ver=version('-release');
        if ~strcmp(ver,'2011b')
            if str2num(ver(1:end-1))<=2011 %#ok<ST2NM>
                xcorrByseg=true;
            end
        end
    end
    if xcorrByseg
        trace1=trace;
        difs=[];
        while length(trace1)>length(HB)*3 %2^20
            if length(trace1)>2^20
                trace2=trace1(1:2^20);
                trace1=trace1(2^20+1:end);
            else
                trace2=trace1;
                trace1=0;
            end
            xcCurrent=xcorr(trace2,HB);
            xcCurrent=xcCurrent(find(xcCurrent,1):end);
            [~,ipxcCur]=findPeaks(xcCurrent,1.5,sRate*period*0.6);
            difs=[difs,diff(ipxcCur)]; %#ok<AGROW>
        end
        period4=median(difs(difs/sRate<maxPeriod))/sRate;
    else
        HBxc1=xcorr(trace,HB);
        [~,ipxc]=findPeaks(HBxc1,1.5,sRate*period*0.6); % index peak xcorr
        period4=median(diff(ipxc))/sRate;
    end
else
    %    warning('could not find cross correlation within extended template, guessing period')
    period4=period;
end
temp=HB(sampBefore-round(sRate*beforeHBs*period4):sampBefore+round(sRate*afterHBs*period4));
edgeRepressor=ones(size(temp));

reducVec=0:1/repressTime:1;
reducVec=reducVec(1:end-1);
edgeRepressor(1:length(reducVec))=reducVec;
edgeRepressor(end-length(reducVec)+1:end)=fliplr(reducVec);
tempe=temp-median(temp);
tempe=tempe.*edgeRepressor;
time=1/sRate:1/sRate:length(temp)/sRate;
time=time-(1-afterHBs)*period4;
if figs
    figure;
    plot(time,tempe,'g')
    title('template HB')
end
tempe=double(tempe);

function MCG=makeMCG(temp,maxTemp,Rlims,Ipeaks,amp,lengt)
MCG=zeros(1,lengt);
HBol=[]; %overlap
olCount=0;
%[~,maxTemp]=max(temp(1:round(length(temp/2))));
for HBi=1:length(Ipeaks);
    s0=Ipeaks(HBi)-maxTemp+1;
    s1=Ipeaks(HBi)+length(temp)-maxTemp;
    if sum(MCG(s0:s1))>0
        overlap=find(MCG(s0:s1),1,'last');
        if overlap>0.2*length(temp)
            endPrev=round(0.2*length(temp));
        else
            endPrev=overlap;
        end
        %         sampDif=round(sRate/50);
        reducVec=1:-1/endPrev:1/endPrev;
        reducVecLR=fliplr(reducVec);
        reducVec(end+1:length(temp))=0;
        reducVecLR(end+1:length(temp))=1;
        MCG(s0:s1)=MCG(s0:s1).*reducVec+temp.*reducVecLR;
        olCount=olCount+1;
        HBol(olCount)=HBi;
    else
        MCG(s0:s1)=MCG(s0:s1)+temp;
    end
    MCG(s0+Rlims(1)-1:s0+Rlims(2)-1)=temp(Rlims(1):Rlims(2))*amp(HBi);
end
% if ~isempty(HBol)
% diary('HBlog.txt')
% disp(['overlapping heartbeats at ',num2str(Ipeaks(HBol/sRate)),'s'])
% diary off
% end
function tempe=HBbyChan(trace,peakIndex,sampBefore,sampAfter,repressTime)
HB=zeros(size(trace,1),sampBefore+1+sampAfter);
% average HBs
for epochi=1:length(peakIndex)
    HB=HB+trace(:,peakIndex(epochi)-sampBefore:peakIndex(epochi)+sampAfter);
end
HB=HB/epochi;
% reduce edges to zero
edgeRepressor=ones(1,size(HB,2));
reducVec=0:1/repressTime:1;
reducVec=reducVec(1:end-1);
edgeRepressor(1:length(reducVec))=reducVec;
edgeRepressor(end-length(reducVec)+1:end)=fliplr(reducVec);
tempe=HB-repmat(mean(HB(:,[1:repressTime,end-repressTime:end]),2),1,size(HB,2));
tempe=tempe.*repmat(edgeRepressor,size(HB,1),1);
function MCG=makeMCGbyCh(temp,maxTemp,Rlims,Ipeaks,amp,lengt)
MCG=zeros(size(temp,1),lengt);
for HBi=1:length(Ipeaks);
    s0=Ipeaks(HBi)-maxTemp+1;
    s1=Ipeaks(HBi)+length(temp)-maxTemp;
    if sum(MCG(1,s0:s1))>0
        overlap=find(MCG(1,s0:s1),1,'last');
        if overlap>maxTemp/2; % 0.2*size(temp,2)
            endPrev=round(maxTemp/2); % round(0.2*size(temp,2));
        else
            endPrev=overlap;
        end
        %         sampDif=round(sRate/50);
        reducVec=1:-1/endPrev:1/endPrev;
        reducVecLR=fliplr(reducVec);
        reducVec(end+1:size(temp,2))=0;
        reducVecLR(end+1:size(temp,2))=1;
        reducVec=repmat(reducVec,size(temp,1),1);
        reducVecLR=repmat(reducVecLR,size(temp,1),1);
        MCG(:,s0:s1)=MCG(:,s0:s1).*reducVec+temp.*reducVecLR;
    else
        MCG(:,s0:s1)=MCG(:,s0:s1)+temp;
    end
    MCG(:,s0+Rlims(1)-1:s0+Rlims(2)-1)=temp(:,Rlims(1):Rlims(2))*amp(HBi);
end
function xcr=XCORR(x,y)
xcorrByseg=false;
if length(x)>=2^20 % test if version is later than 1011b
    ver=version('-release');
    if ~strcmp(ver,'2011b')
        if str2num(ver(1:end-1))<=2011
            xcorrByseg=true;
        end
    end
end
[~,Rsamp]=max(y);
if xcorrByseg
    trace1=x;
    xcr=[];
    
    while length(trace1)>length(y)%2^20
        if length(trace1)>2^20
            trace2=trace1(1:2^20);
            trace1=trace1(2^20+1:end);
        else
            trace2=trace1;
            trace1=0;
        end
        xcrPad=zeros(size(trace2));
        [xcCurrent,lags]=xcorr(trace2,y);
        xcCurrent=xcCurrent(lags>=0);
        xcrPad(Rsamp:end)=xcCurrent(1:end-Rsamp+1);
        xcr=[xcr,xcrPad];
        % FIXME xcorr (fftfilt) won't take it for more than 2^20
    end
    
else
    [xcr,lags]=xcorr(x,y);
    xcr=xcr(lags>=0);
    xcrPad=zeros(size(x));
    xcrPad(Rsamp:end)=xcr(1:end-Rsamp+1);
    xcr=xcrPad; % sorry for switching variables
end
function [p,Rlims]=assessAmp(templateHB,maxi,Ipeaks2in,meanMEG)
%[~,maxi]=max(templateHB(1:round(length(templateHB/2))));
bef=find(fliplr(templateHB(1:maxi))<=0,1)-1;
aft=find(templateHB(maxi:end)<=0,1)-1;
Rlims=[maxi-bef,maxi+aft]; % check where R pulls the template above zero
for HBi=1:length(Ipeaks2in);
    s0=Ipeaks2in(HBi)-bef;
    s1=Ipeaks2in(HBi)+aft;
    x=templateHB(Rlims(1):Rlims(2));
    y=meanMEG(s0:s1);
    scalef=-round(log10(max([x,y]))); % scaling factor for polyfit not to complain
    x=x*10^scalef;
    y=y*10^scalef;
    pScaled=polyfit(x,y,1);
    pScaled(2)=pScaled(2).*10^-scalef;
    p(HBi,1:2)=pScaled; %#ok<AGROW>
end
function variable=default(field,value,cfg)
if isfield(cfg,field)
    eval(['variable=cfg.',field,';'])
else
    variable=value;
end
function vecN=standard(vec) % normalize data for display
vecN=vec/median(abs(vec));
