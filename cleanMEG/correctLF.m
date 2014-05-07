function [cleanData,whereUp,noiseSamp]=correctLF(data,sRate,chanLF,method,Lfreq,jobs,hpFreq,noiseThr)
%   -  data is MEG or EEG data, rows for channels
%   -  sRate is sampling rate
%   -  chanLF is a channel containing a lot of 50 or 60Hz, can be a MEG
% reference channel or leave empty to detect high LF chan automatically.
% you can also give 'time', the consequence is that LF will not be searched
% on any channel, but will be guessed according to time only. requires Lfreq and works with ADAPTIVE method.
%   -  method is 'GLOBAL' or 'ADAPTIVE' (default), see cleanLineF for more on that.
%   -  Lfreq is the noise frequency, 50 or 60. by default the Lfreq is
%   detected automatically
%   - jobs is for parallel processing with parfor, a number of CPU to use.
%   default uses all, take care.
%   - hpFreq (0.1) means that you want to highpass filter the data before
%   cleaning.
%   - noiseThr is how many std(mean(abs(data of one cycle))) to consider as
%   noise. default is 5 SD.

% 4D users can run it as cleanData=LFcleanNoCue;
%
% Yuval Harpaz Jan 2014

%% set defaults
if ~exist ('method','var')
    method=[];
end
if isempty(method)
    method='ADAPTIVE';
end
if ~ischar(method)
    Ncycles=method;
    method='ADAPTIVE';
else
    Ncycles=[];
end
if ~exist ('chanLF','var')
    chanLF=[];
end
if ~exist('Lfreq','var')
    Lfreq=[];
end
if ~exist('jobs','var')
    jobs=[];
end
if ~exist('hpFreq','var')
    hpFreq=[];
end
if ~exist('noiseThr','var')
    noiseThr=[];
else
    if ~strcmp(method,'ADAPTIVE')
        warning('currently noiseThr only works for adaptive')
    end
end
if isempty(noiseThr)
    noiseThr=5;
end
%% try to load data file and check if 4D-neuroimaging data

if ~exist('data','var')
    data=[];
    sRate=[];
end
if ischar(data)
    if strcmp(data(end-3:end),'.mat') % read matrix from file 'data.mat'
        PWD=pwd;
        display(['loading ',PWD,'/',data,]);
        data=load(['./',data]);
        dataField=fieldnames(data);
        eval(['data=data.',dataField{1,1},';']);
    else % read 4D data from file name specified in 'data'
        cloc=strfind(data,'c');
        comaloc=strfind(data,',');
        if ~isempty(cloc) && ~isempty(comaloc)
            if comaloc(end)>cloc(1)
                var4DfileName=data;
            end
        end
    end
end

if isempty(data) || exist('var4DfileName','var');
    if ~exist('var4DfileName','var');
        try
            var4DfileName=ls('hb_c,*');
        catch %#ok<CTCH>
            var4DfileName=ls('c,*');
        end
        var4DfileName=['./',var4DfileName(1:end-1)];
    end
    var4Dp=pdf4D(var4DfileName);
    sRate=double(get(var4Dp,'dr'));
    var4Dhdr = get(var4Dp, 'header');
    var4DnSamp=var4Dhdr.epoch_data{1,1}.pts_in_epoch;
    var4Dchi = channel_index(var4Dp, 'meg', 'name');
    display(['reading ',var4DfileName]);
    data = read_data_block(var4Dp,[1 var4DnSamp],var4Dchi);%,var4Dchi);
    if isempty(chanLF)
        var4Dlabels={'MCxaA','MCyaA','MCzaA','MLxaA','MLyaA','MLzaA','MRxaA','MRyaA','MRzaA'};
        var4DchiRef = channel_index(var4Dp,var4Dlabels, 'name');
        var4Dref = read_data_block(var4Dp,[1 var4DnSamp],var4DchiRef);%,var4Dchi);
        
        [FourRef,Fref]=fftBasic(var4Dref,round(sRate));
        if isempty(Lfreq)
            Lfreq=findLfreq(FourRef,Fref);
        end
        if isempty(chanLF)
            maxChani=findchanLF(FourRef,Fref,Lfreq);
        end
        chanLF=var4Dref(maxChani,:);
        display(['selected ',var4Dlabels{maxChani},' as Line Freq cue'])
    end
    clear var4D*
end
if isempty(jobs) || size(data,1)==1
    par=false;
else
    par=true;
end
%% find chans with no obvious problems and compute meanPSD

testSamp=min([round(sRate) size(data,2)]);
for chani=1:size(data,1)
    good(chani)=true; %#ok<AGROW>
    if isequal(data(chani,1:testSamp),int16(data(chani,1:testSamp))) || length(unique(data(chani,1:testSamp)))<20
        good(chani)=false; %#ok<AGROW>
    end
end
good=find(good);
%% if requested make hp filter
if ~isempty(hpFreq);
    display('filtering')
    HighPassSpecObj=fdesign.highpass('Fst,Fp,Ast,Ap',0.001,hpFreq,60,1,sRate);%
    HighPassFilt=design(HighPassSpecObj ,'butter');
    % BL correction for 1st sec window to avoid ripples in the beginning
    for chani=1:size(data,1)
        data(chani,:)=data(chani,:)-mean(data(chani,1:round(sRate)));
    end
    data = myFilt(data,HighPassFilt);
end
%% test 50Hz or 60Hz
display('computing fft')
% check size not to overload
dSize=size(data,1)*size(data,2);
if dSize>500000000
    fftLength=round(500000000/size(data,1));
    display(['lots of data, displaying fft only for about ',num2str(round(fftLength/sRate/60)),'min'])
else
    fftLength=size(data,2);
end
[Four,F]=fftBasic(data(good,1:fftLength),round(sRate));
if isempty(Lfreq);
    display('filtering and searching for line frequency')
    [Lfreq,meanPSD]=findLfreq(Four,F);
    disp(['line frequency is ',num2str(Lfreq),'Hz'])
else
    [~, i125] = min(abs(F-125)); % index for 125Hz
    [~, i145] = min(abs(F-145));
    scale=mean(abs(Four(:,i125:i145)),2);
    FourScaled=zeros(size(Four));
    for chani=1:size(Four,1)
        FourScaled(chani,:)=abs(Four(chani,:))/scale(chani);
    end
    meanPSD=mean(FourScaled,1);
end
display('done fft, preparing LF time indices')
%% find LF cycles on a filtered channel
lookForLF=true;

if  isempty(chanLF);
    %[Four,F]=fftBasic(data(good,:),round(sRate));
    
    maxChani=findchanLF(Four,F,Lfreq);
    %     [~, i125] = min(abs(F-125)); % index for 125Hz
    %     [~, i145] = min(abs(F-145));
    %     [~, iLF] = min(abs(F-Lfreq));
    %     scale=mean(abs(Four(:,i125:i145)),2);
    %     [~,maxChani]=max(abs(Four(:,iLF))./scale); %round(freq)==Lfreq
    chanLF=data(good(maxChani),:);
    display(['selected chan index ',num2str(good(maxChani)),' as Line Freq cue'])
elseif ischar(chanLF)
    if strcmp(chanLF,'time')
        % cycleSamp=round(sRate/Lfreq)
        start=double(1):sRate/Lfreq:double(size(data,2));
        whereUp=round(start);
        lookForLF=false;
    elseif strcmp(chanLF(end-3:end),'.mat')
        chanLF=load(['./',data]);
        dataField=fieldnames(chanLF);
        eval(['chanLF=chanLF.',dataField{1,1},';']);
    else
        error('could not figure out what to do with chanLF string')
    end
elseif length(chanLF)==1;
    chanLF=data(chanLF,:);
end
if lookForLF
    if length(unique(chanLF))==2 % zero crossing is coded in chanLF (trigger channel)
        trigShift=zeros(size(chanLF));
        trigShift(2:end)=chanLF(1:end-1);
        trigShift=uint16(trigShift);
        whereUp=find((uint16(chanLF)-trigShift)>0);
        if max(diff(whereUp))>sRate/Lfreq+3 || min(diff(whereUp))<sRate/Lfreq-3
            warning('some triggers are at irregular times')
        end
    else
        if ~exist('BPfilt','var')
            BPobj=fdesign.bandpass(...
                'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
                Lfreq-10,Lfreq-5,Lfreq*2+5,Lfreq*2+10,60,1,60,sRate);
            BPfilt=design(BPobj ,'butter');
        end
        chanLFf=myFilt(chanLF-mean(chanLF(1:round(sRate))),BPfilt);
        chanLFShift=[0,chanLFf(1:end-1)];
        logsum=chanLFf>0;
        logsum2=chanLFShift<0;
        whereUp=logsum+logsum2==2;
        whereUp=find(whereUp);
    end
end
%% check that whereUp is OK. if not, find cycles by template matching
cycInterval=1000*mean(diff(whereUp))/sRate;
if ~round(cycInterval)==1000/Lfreq
    warning('whereUp has wrong frequency')
end
badCue=find(diff(whereUp)>1.2*sRate/Lfreq); %#ok<EFIND>
if ~isempty(badCue)
    % try matching a template
    hpObj=fdesign.highpass('Fst,Fp,Ast,Ap',0.1,round(Lfreq*0.67),60,1,sRate);%
    Filt=design(hpObj ,'butter');
    chanLFhp = myFilt(chanLF,Filt);
    cycCount=0;
    temp=zeros(1,round(sRate/Lfreq));
    for upi=round(1:sRate/Lfreq:sRate*5)
        cycCount=cycCount+1;
        temp(cycCount,:)=chanLFhp(upi:upi+round(sRate/Lfreq)-1);
    end
    temp=mean(temp,1);
    [snr,signal]=match_temp(chanLFhp,temp,1);
    [~, IpeaksPos]=findPeaks(signal,0,round(0.75*length(temp))); % no peaks closer than 60% of period
    [~, IpeaksNeg]=findPeaks(-signal,0,round(0.75*length(temp)));
    
    %     [~, IpeaksPos]=findPeaks(snr,0,round(0.75*length(temp))); % no peaks closer than 60% of period
    %     [~, IpeaksNeg]=findPeaks(-snr,0,round(0.75*length(temp))); % no peaks closer than 60% of period
    if std(diff(IpeaksPos))<std(diff(IpeaksNeg))
        Ipeaks=IpeaksPos;
    else
        Ipeaks=IpeaksNeg;
    end
    if std(diff(Ipeaks))>0.2*sRate/Lfreq
        warning('std of cue interval is large')
    end
    
    % whereUp=setdiff(Ipeaks,Ipeaks(badCue)); % to ignore bad cues
    whereUp=Ipeaks;
    
end
cycInterval=1000*mean(diff(whereUp))/sRate;
if ~round(cycInterval)==1000/Lfreq
    error('whereUp has wrong frequency')
end




if whereUp(1)>(ceil(sRate/Lfreq)+1)
    diary LFlog.txt
    warning('first LF index is away from the beginning, guessing first cycles')
    diary off
    firstInd=fliplr(round(double(whereUp(1)):-mean(diff(whereUp)):1));
    whereUp=[firstInd(1:end-1),whereUp];
end
if (size(data,2)-whereUp(end))>(ceil(sRate/Lfreq)+1)
    diary LFlog.txt
    warning('last LF index is away from the end, guessing last cycles')
    diary off
    lastInd=round(double(whereUp(end)):mean(diff(whereUp)):size(data,2));
    whereUp=[whereUp,lastInd(2:end)];
end
badCue=find(diff(whereUp)>ceil(1.2*sRate/Lfreq));
if ~isempty(badCue)
    diary LFlog.txt
    disp(['irregular interval between cues at (seconds): ',num2str(whereUp(badCue)/sRate)])
    diary off
end

%% prepare parallel processing
if par
    closeLabs=false;
    try
        nCPU=matlabpool('size'); % 0 if no matlabpool yet;
    catch
        disp('could not find "matlabpool". Setting "nCPU" to 1')
        nCPU=1;
    end
    if nCPU==0
        closeLabs=true;
        if isempty(jobs)
            try
                matlabpool;
                jobs=matlabpool('size');
            catch
                jobs=1;
                warning('cannot get matlabpool to work, define parallel for your matlab to work faster')
            end
        else
            try
                matlabpool ('open', jobs);
            catch
                jobs=1;
                warning('cannot get matlabpool to work, define parallel for your matlab to work faster')
            end
        end
    else
        if isempty(jobs)
            try
                jobs=matlabpool('size');
            catch
                jobs=1;
                warning('cannot get matlabpool to work, define parallel for your matlab to work faster')
            end
        elseif nCPU~=jobs
            try
                matlabpool close
                %matlabpool('local',jobs);
                matlabpool ('open', jobs)
            catch
                jobs=1;
                warning('cannot get matlabpool to work, define parallel for your matlab to work faster')
            end
        end
    end
end
%% cleaning the data
% estimate time
if size(data,1)==1
    [cleanData,~,nS]=cleanLineF(data, whereUp, [], upper(method),[],Ncycles,noiseThr);
    noiseSamp{1,1}=nS;
    clear data
else
    tic
    cleanLineF(data(good(1),:), whereUp, [], upper(method),[],Ncycles,noiseThr);
    time1st=toc;
    if par
        timeEstimate=num2str(ceil(time1st*(length(good)-1)/60/jobs*2));
        display(['cleaning channels one by one, wait about ',timeEstimate,'min'])
    else
        timeEstimate=num2str(ceil(time1st*(length(good)-1)/60/1*2));
        display(['cleaning channels one by one, wait about ',timeEstimate,'min'])
    end
    % do the cleaning
    cleanData=data(good,:);
    nSamp=cell(size(cleanData,1),1);
    if par
        parfor chani=1:length(good)
            [cleanData(chani,:),~,nSamp{chani,1}]=cleanLineF(cleanData(chani,:), whereUp, [], upper(method),[],Ncycles,noiseThr);
        end
    else
        for chani=1:length(good)
            [cleanData(chani,:),~,nSamp{chani,1}]=cleanLineF(cleanData(chani,:), whereUp, [], upper(method),[],Ncycles,noiseThr);
        end
    end
    data(good,:)=cleanData;
    clear cleanData
    cleanData=data;% sorry about this mess, the parfor made me do this
    clear data
    noiseSamp=cell(size(cleanData,1),1);
    for chani=1:length(good)
        noiseSamp{good(chani),1}=nSamp{chani,1};
    end
    time2nd=toc;
    mins=time2nd/60;
    if floor(mins)>str2num(timeEstimate)
        well='well, ';
    else
        well='';
    end
    secs=60*(mins-floor(mins));
    display([well,'it took ',num2str(floor(mins)),'min and ',num2str(round(secs)),'sec']);
end
disp('plotting')
[Four,F]=fftBasic(cleanData(good,1:fftLength),round(sRate));
[~, i125] = min(abs(F-125)); % index for 125Hz
[~, i145] = min(abs(F-145));
scale=mean(abs(Four(:,i125:i145)),2);
for chani=1:size(Four,1)
    Four(chani,:)=abs(Four(chani,:))/scale(chani);
end
meanPSDclean=mean(Four,1);
figure;
plot(F,meanPSD,'r');
hold on
plot(F,meanPSDclean,'g')
legend('original','clean')
title('PSD after rescaling and averaging channels')
if par
    if closeLabs
        try
            matlabpool close
        end
    end
end
display('done cleaning LF')
function [Lfreq,meanPSD]=findLfreq(fourier,freq)
% [fourier,freq]=fftBasic(data,round(sRate));
[~, i125] = min(abs(freq-125)); % index for 125Hz
[~, i145] = min(abs(freq-145));
scale=mean(abs(fourier(:,i125:i145)),2);
for chani=1:size(fourier,1)
    fourier(chani,:)=abs(fourier(chani,:))/scale(chani);
end
meanPSD=mean(fourier,1); % power spectrum, averaged over channels
[~, i50] = min(abs(freq-50));
[~, i60] = min(abs(freq-60));

snr50=2*meanPSD(i50)/(meanPSD(i50-2)+meanPSD(i50+2));
snr60=2*meanPSD(i60)/(meanPSD(i60-2)+meanPSD(i60+2));
if meanPSD(i50)>meanPSD(i60) && snr50>2
    Lfreq=50;
elseif meanPSD(i60)>meanPSD(i50) && snr60>2
    Lfreq=60;
else
    plot(freq,meanPSD)
    title('Power Spectrum averaged over channels')
    error('cannot makeup my mind if thhere is 50 or 60 Hz artifact')
end
function maxChani=findchanLF(fourier,freq,Lfreq)
[~, i125] = min(abs(freq-125)); % index for 125Hz
[~, i145] = min(abs(freq-145));
[~, iLF] = min(abs(freq-Lfreq));
scale=mean(abs(fourier(:,i125:i145)),2);
[~,maxChani]=max(abs(fourier(:,iLF))./scale); %round(freq)==Lfreq



