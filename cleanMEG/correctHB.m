function [cleanData,temp2e,period2,ECG,Rtopo]=correctHB(data,sRate,figOptions)%rawData,sRate,ecg)
% data is a matrix with rows for channels, raw data, not filtered.
% sRate is sampling rate
% figs=false;
% if you want topoplot of HB (first and second sweeps) you have to have:
% figOptions.label=var4Dlabel;
% figOptions.layout='4D248.lay';
% 4D users can run the function from the folder with 'c,*' file with no
% input arguments:
% cleanData=correctHB;
% or like this cleanData=correctHB([],[],1); to get the figures.
% if you don't specify figure options you still get one before / after figure.
% added by Dr. Yuval Harpaz to Prof. Abeles' work

% Issues
%  - no fail warning for data with no HB artifact
%  - overlap between two heartbeats may have some nosie from summing two
% templates
%  - designed for magnetometers, needs adjustment for mag + grad data, like
%  finding peaks using magnetometers and removing averaged template for
%  each channel (including grad)
%
% it works, try it!

%% default variables and parameters
if ~exist('figOptions','var')
    figOptions=[];
end
if isempty(figOptions)
    figs=false;
else
    figs=true;
end
minPeriod=0.45;
maxPeriod=1.5;
%% checking defaults for 4D data
% to use with data=[] and sRate=[];
if ~exist('data','var')
    data=[];
    sRate=[];
end
if isempty(data);
    var4DfileName=ls('c,*');
    var4DfileName=var4DfileName(1:end-1);
    var4Dp=pdf4D(var4DfileName);
    sRate=double(get(var4Dp,'dr'));
    var4Dhdr = get(var4Dp, 'header');
    var4DnSamp=var4Dhdr.epoch_data{1,1}.pts_in_epoch;
    var4Dchi = channel_index(var4Dp, 'meg', 'name');
    data = read_data_block(var4Dp,[1 var4DnSamp],var4Dchi);
    
    if figs
        var4Dlabel=channel_label(var4Dp,var4Dchi)';
        figOptions.label=var4Dlabel;
        figOptions.layout='4D248.lay';
    end
    clear var4D*
end

%% look for a noisy segment and noisy channels
zThr=20;
% find bad channels, has to be noisy for 3 of the first 3 seconds

badc=zeros(size(data,1),1); % bad channels
stdMEG=std(data(:,1:round(sRate))');
badc=stdMEG>zThr;
stdMEG=std(data(:,round(sRate):round(2*sRate))');
badc=badc+stdMEG>zThr;
stdMEG=std(data(:,round(2*sRate):round(3*sRate))');
badc=badc+stdMEG>zThr;
badc=find(badc==3);
if ~isempty(badc)
    data(badc,:)=0;
end

% find jump or other huge artifact
zThr=10;
meanMEG=mean(data);
meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
zMEG=(meanMEGdt-mean(meanMEGdt))./std(meanMEGdt);
jbeg=find(abs((zMEG))>zThr,1);
bads=[]; % bad samples
if ~isempty(jbeg)
    jend=find(abs((zMEG))>zThr,1,'last');
    jend2=find(abs(zMEG(jend:end))>1,1,'last')+jend-1;
    bads=(jbeg-round(sRate./2)):(jend2+round(sRate*0.5));
    %     if length(bads)>sRate*5
    %         error(['jump? what''','s the noise at ',num2str(jbeg./sRate),'s ? more then 5sec between first and last noisy point, cant handle this, please clean piece by piece']);
    %     else
    data(:,bads)=0;
    data=data-repmat(median(data(:,1:jbeg),2),1,size(data,2));
    diary('HBlog.txt')
    warning(['jump? what''','s the noise at ',num2str(jbeg./sRate),'s? zeroed from ',num2str(bads(1)/sRate),' to ',num2str(bads(end)/sRate)]);
    diary off
    %    end
else
    % baseline correction by removing the median of each channel
    data=data-repmat(median(data,2),1,size(data,2));
end




%% filter mean MEG data



% averaging the MEG channels
meanMEG = mean(data);
meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
% filtering to pass from 5-7Hz to 100-110Hz
BandPassSpecObj=fdesign.bandpass(...
    'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
    5,7,90,110,60,1,60,sRate);
BandPassFilt=design(BandPassSpecObj  ,'butter');
meanMEGf = myFilt(meanMEG,BandPassFilt);
% baseline correction again, just in case
meanMEGf=meanMEGf-median(meanMEGf);

sampBefore=round(sRate*maxPeriod);
%% peak detection on ECG like signal


[peaks, Ipeaks]=findPeaks(meanMEGf,1.5,round(sRate*minPeriod)); % 450ms interval minimum
if figs
    figure;
    plot(meanMEGf)
    hold on
    plot(Ipeaks, peaks,'ro')
    title('peak detection on mean MEG trace, OK if most of them are HB')
end
%% get topography
if figs
    if isfield(figOptions,'layout') && isfield(figOptions,'label')
        topo={};
        topo.avg=median(data(:,Ipeaks),2);
        topo.time=0;
        topo.label=figOptions.label;
        topo.dimord='chan_time';
        cfg=[];
        cfg.layout=figOptions.layout;
        %cfg.interpolation='linear';
        cfg.xlim=[1,1];
        cfg.zlim=[-max(abs(topo.avg)) max(abs(topo.avg))];
        figure;
        ft_topoplotER(cfg,topo);
    else
        warning('no topoplot without labels and layout fields! see figOptions options')
    end
end

topoTrace=median(data(:,Ipeaks),2)'*data;
topoTrace=myFilt(topoTrace,BandPassFilt);
%dataFilt=myFilt(data,BandPassFilt);
topoTrace=topoTrace-median(topoTrace);
topoTraceN=topoTrace./max(topoTrace(1:round(sRate*10)));
meanMEGN=meanMEGf./max(meanMEGf(1:round(sRate*10)));

%% check if topo of every peak is correlated to average topo
r=corr(data(:,Ipeaks),median(data(:,Ipeaks),2));
if figs
    figure;
    plot(topoTraceN)
    hold on
    plot(meanMEGN,'r')
    plot(Ipeaks(r>0.5),r(r>0.5),'g.');
    legend('topoTrace','meanMEG','r data-topo > 0.5')
end
%% average good HB and make a template
IpeaksR=Ipeaks(r>0.5);
IpeaksR=IpeaksR(IpeaksR>sampBefore);
IpeaksR=IpeaksR(IpeaksR<(size(data,2)-sampBefore));

[temp1e,period1]=makeTempHB(mean(data),sRate,IpeaksR,minPeriod/0.6,sampBefore,figs);
%% find xcorr between template and meanMEG
[xcr,lags]=xcorr(meanMEGf,temp1e);
% figure;
% plot(lags,xcr);
xcr=xcr(lags>=0);
lags=lags(lags>=0);
[~,tempMax]=max(temp1e);
%xr(length(t)+1:end)=x(1:end-length(t));

xcrPad=zeros(size(meanMEGf));

xcrPad(tempMax:end)=xcr(1:end-tempMax+1);
% figure;
% plot(meanMEGN);
% hold on
% plot(xcrPad/max(xcrPad),'r');
if figs
    figure;
    plot(topoTraceN)
    hold on
    plot(meanMEGN,'r')
    plot(xcrPad/max(xcrPad),'g');
    legend('topoTrace','meanMEG','temp xcorr')
end
%% second sweep
%% find peaks on xcorr trace
[peaks2, Ipeaks2]=findPeaks(xcrPad,1.5,round(sRate*period1*0.6)); % no peaks closer than 60% of period
if figs
    figure;
    plot(xcrPad)
    hold on
    plot(Ipeaks2, peaks2,'ro')
    title('2nd sweep peak detection, based on template matching')
end
% checking results
periodS=diff(Ipeaks2);% the period in samples for each HB
farI=find(periodS/median(periodS)>1.5);
if median(periodS)/sRate<0.5 || median(periodS)/sRate>1.5
    diary('HBlog.txt')
    warning('interval between HB not resonable, possibly,not enogh HBs detected. look at the figures!')
    diary off
else
    if ~isempty(farI)
        farT=Ipeaks2(farI)/sRate; % this is far time, not what you think
        diary('HBlog.txt')
        disp('sparse heartbeats, you may want to look for missing HB after:');
        disp(num2str(farT));
        diary off
    end
    nearI=find(periodS/median(periodS)<0.5);
    
    if ~isempty(nearI)
        nearT=Ipeaks2(nearI)/sRate; % this is far time, not what you think
        diary('HBlog.txt')
        disp('close heartbeats, you may want to look for false HB detection after:');
        disp(num2str(nearT));
        diary off
    end
end

% ignore edges
Ipeaks2in=Ipeaks2(Ipeaks2>sampBefore);
Ipeaks2in=Ipeaks2in(Ipeaks2in<(size(data,2)-sampBefore));
% make ecg trace for meanMEG
[temp2e,period2]=makeTempHB(meanMEG,sRate,Ipeaks2in,period1,sampBefore,figs);
% FIXME reject bad SNR HB from averaging



%% test R amplitude
% meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
[~,maxi]=max(temp2e(1:round(length(temp2e/2))));
bef=find(flipLR(temp2e(1:maxi))<0,1)-1;
aft=find(temp2e(maxi:end)<0,1)-1;
Rlims=[maxi-bef,maxi+aft]; % check where R pulls the template above zero

% ampMM=meanMEGdt(Ipeaks2in)./temp2e(maxi);
% mmSm=smooth(meanMEGdt,10);
% tSm=smooth(temp2e,10);
% ampMMsm=mmSm(Ipeaks2in)./tSm(maxi);

for HBi=1:length(Ipeaks2in);
    s0=Ipeaks2in(HBi)-bef;
    s1=Ipeaks2in(HBi)+aft;
    p(HBi,1:2)=polyfit(temp2e(Rlims(1):Rlims(2)),meanMEGdt(s0:s1),1);
end
negp=find(p(:,1)<0);
if ~isempty(negp)
    p(negp,1:2)=0;
end
ampMMfit=p(:,1)+(p(:,2)./temp2e(maxi));

ECG=makeECG(temp2e,Rlims,Ipeaks2in,ampMMfit,length(meanMEG));
% figure;
% plot(meanMEGdt,'r')
% hold on
% plot(ECG)
% legend('detrended meanMEG','line fit amp ECG')


%% remove ecg from each chan


%make temp per chan
sBef=maxi-1;
sAft=length(temp2e)-maxi;
HBtemp=HBbyChan(data,sRate,Ipeaks2in,period2,sBef,sAft);
% check SNR per chan
sm50=round(sRate*0.05);
s0=maxi-sm50*3;
s1=maxi-sm50;
s2=maxi+sm50;
%HBsnr=HBtemp(:,maxi)./max(abs(HBtemp(:,1:round(sBef/2))),[],2);
n=std(HBtemp(:,s0:s1)');
s=std(HBtemp(:,s1:s2)');
snr=s./n;
badSNR=find(snr<=2);
if figs
    figure;plot(HBtemp','k');hold on;
    %plot(maxi,HBtemp(find(HBsnr>1.1),maxi),'g.')
    [~,minI]=min(snr);
    plot(HBtemp(minI,:),'r');
    if isempty(badSNR)
        title('red trace is worst SNR channel')
    else
        plot(maxi,HBtemp(badSNR,maxi),'r.')
        title('red dots mark channels with SNR < 2')
    end
end

if ~isempty(badSNR)
    if length(find(snr<=2))>size(data,1)/2
        error('too many channels have poor SNR')
    end
    HBtemp(badSNR,:)=0;
    diary('HBlog.txt')
    disp(['not including bad SNR channels, index no. ',num2str(badSNR)])
    diary off
end
% FIXME check poor SNR chans on chan temp
ECGall=makeECGbyCh(HBtemp,Rlims,Ipeaks2in,ampMMfit,length(meanMEG),maxi);

cleanData=data-ECGall;
figure;
plot(ECG,'k')
hold on
plot(meanMEG,'r')
plot(mean(cleanData),'g')
legend('ECG from template', 'mean MEG','mean clean MEG')
Rtopo=HBtemp(:,maxi);
if figs
    topo.avg=Rtopo;%median(dataFilt.trial{1,1}(:,Ipeaks),2);
    cfg.zlim=[-max(abs(Rtopo)) max(abs(Rtopo))];
    figure;
    ft_topoplotER(cfg,topo);
    title ('TOPOGRAPHY OF R')
end

% FIXME check missing / extra HB

function [tempe,period2]=makeTempHB(trace,sRate,peakIndex,period,sampBefore,figs)
betweenHBs=0.7; % after the T wave, before next qrs, 0.7 of period
HB=zeros(size(trace,1),sampBefore*2+1);
for epochi=1:length(peakIndex)
    HB=HB+trace(:,peakIndex(epochi)-sampBefore:peakIndex(epochi)+sampBefore);
end
HB=HB/epochi;
period2=[];

HBxc=xcorr(HB);
% FIXME admit defeat for yoni
[~,ipxc]=findPeaks(HBxc,1.5,sRate*period*0.6); % index peak xcorr
if length(ipxc)>1
    nextPeak=ceil(length(ipxc)/2+0.5);
    period2=(ipxc(nextPeak)-ipxc(nextPeak-1))/sRate;
else
    HBxc1=xcorr(trace,HB);
    [~,ipxc]=findPeaks(HBxc1,1.5,sRate*period*0.6); % index peak xcorr
    
    period2=median(diff(ipxc))/sRate;
end

display(['HB period (2nd sweep) is ',num2str(period2),'s']);

temp=HB(sampBefore-round(sRate*(1-betweenHBs)*period2):sampBefore+round(sRate*betweenHBs*period2));
edgeRepressor=ones(size(temp));
ms20=round(sRate/50);
reducVec=[0:1/ms20:1];
reducVec=reducVec(1:end-1);
edgeRepressor(1:length(reducVec))=reducVec;
edgeRepressor(end-length(reducVec)+1:end)=flipLR(reducVec);
tempe=temp-median(temp);
tempe=tempe.*edgeRepressor;
if figs
    figure;
    plot(tempe,'g')
    title('template HB')
end
tempe=double(tempe);

function ECG=makeECG(temp,Rlims,Ipeaks,amp,lengt)
ECG=zeros(1,lengt);
[~,maxTemp]=max(temp(1:round(length(temp/2))));
% bef=find(flipLR(temp(1:maxTemp))<0,1)-1;
% aft=find(temp(maxTemp:end)<0,1)-1;
% Rlims=[maxTemp-bef,maxTemp+aft];

for HBi=1:length(Ipeaks);
    s0=Ipeaks(HBi)-maxTemp+1;
    s1=Ipeaks(HBi)+length(temp)-maxTemp;
    ECG(s0:s1)=temp;
    ECG(s0+Rlims(1)-1:s0+Rlims(2)-1)=temp(Rlims(1):Rlims(2))*amp(HBi);
end
function tempe=HBbyChan(trace,sRate,peakIndex,period,sampBefore,sampAfter)
betweenHBs=0.7; % after the T wave, before next qrs, 0.7 of period
HB=zeros(size(trace,1),sampBefore+1+sampAfter);
% average HBs
for epochi=1:length(peakIndex)
    HB=HB+trace(:,peakIndex(epochi)-sampBefore:peakIndex(epochi)+sampAfter);
end
HB=HB/epochi;
% reduce edges to vero
edgeRepressor=ones(1,size(HB,2));
ms20=round(sRate/50);
reducVec=[0:1/ms20:1];
reducVec=reducVec(1:end-1);
edgeRepressor(1:length(reducVec))=reducVec;
edgeRepressor(end-length(reducVec)+1:end)=flipLR(reducVec);
tempe=HB-repmat(mean(HB(:,[1:ms20,end-ms20:end]),2),1,size(HB,2));
tempe=tempe.*repmat(edgeRepressor,size(HB,1),1);
function ECG=makeECGbyCh(temp,Rlims,Ipeaks,amp,lengt,maxTemp)
ECG=zeros(size(temp,1),lengt);
for HBi=1:length(Ipeaks);
    s0=Ipeaks(HBi)-maxTemp+1;
    s1=Ipeaks(HBi)+length(temp)-maxTemp;
    ECG(:,s0:s1)=temp;
    ECG(:,s0+Rlims(1)-1:s0+Rlims(2)-1)=temp(:,Rlims(1):Rlims(2))*amp(HBi);
end

