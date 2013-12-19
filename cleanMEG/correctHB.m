function [cleanData,temp2e,period2,MCG,Rtopo]=correctHB(data,sRate,figOptions,ECG,snrThr)

% - data is a matrix with rows for channels, raw data, not filtered. it can
% also be a filename.mat, directing to data matrix file, or a 4D filename such
% as 'c1,rfhp0.1Hz'.
% - sRate is sampling rate
% - figOptions=false;
% if you want fieldtrip topoplot of HB (first and second sweeps) you have to have:
% figOptions.label=N by 1 cell arraay with channel names of data
% figOptions.layout='4D248.lay' for 4D users. I recommend
% 'neuromag306mag.lay' for neuromag users even if data includes also grads.
% - ECG can be ECG (useful for ctf users) or a mean of subset of MEG channels where the HB is
% visible. neuromag users can put there the mean of the magnetometers to
% clean both magnetometers and gradiometers included in data.
% - nsrThr is the threshold that tell which channels are cleaned and which
% remain as are. use 0 to clean all, default is 1.5
%
% 4D users can run the function from the folder with 'c,*' file with no
% input arguments:
% cleanData=correctHB;
% or like this cleanData=correctHB([],[],1); to get the figures.
% if you don't specify figure options you still get one before / after figure.
% added by Dr. Yuval Harpaz to Prof. Abeles' work

% Issues
% - amplitude estimate per HB needs be better testing
% - only R amplitude is corrected, may consider to change q s and t waves.
% - memory problem, I should fix it to make topo and template based on the
% beginning of the data to clean the rest with it as well, piece by piece.
% - allow using a premade template, good for cleaning data in 2 pieces
% - 
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
maxPeriod=2; % M sometimes had 1.8s period
if ~exist('ECG','var')
    ECG=[];
end
if ~exist('snrThr','var')
    snrThr=[];
end
if isempty(snrThr)
    snrThr=1.5;
end
%% checking defaults for 4D data
% to use with data=[] and sRate=[];
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
            var4DfileName=ls('xc,lf_c,*');
        catch
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
    data = read_data_block(var4Dp,[1 var4DnSamp],var4Dchi);
    %data=double(data);
    if figs
        var4Dlabel=channel_label(var4Dp,var4Dchi)';
        figOptions.label=var4Dlabel;
        figOptions.layout='4D248.lay';
    end
    clear var4D*
end
%% filter mean MEG (or ECG) data
if ~isempty(ECG)
    meanMEG=ECG;
    %meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
else
    meanMEG=double(mean(data));
end
% filtering to pass from 5-7Hz to 100-110Hz
BandPassSpecObj=fdesign.bandpass(...
    'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
    5,7,90,110,60,1,60,sRate);
BandPassFilt=design(BandPassSpecObj ,'butter');
meanMEGf = myFilt(meanMEG,BandPassFilt);
% baseline correction again, just in case
meanMEGf=meanMEGf-median(meanMEGf);
sampBefore=round(sRate*maxPeriod);
%% look for a noisy segment and noisy channels
zThr=20;
% find bad channels, has to be noisy for 3 of the first 3 seconds
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
zThr=15;
%meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
zMEG=(meanMEGf-mean(meanMEGf))./std(meanMEGf);
jbeg=find(abs((zMEG))>zThr,1);
bads=[]; %#ok<NASGU> % bad samples
if ~isempty(jbeg)
    jend=find(abs((zMEG))>zThr,1,'last');
    jend2=find(abs(zMEG(jend:end))>1,1,'last')+jend-1;
    bads=(jbeg-round(sRate./2)):(jend2+round(sRate*0.5));
    data(:,bads)=0;
    if length(data)<2^19
        data=data-repmat(median(data(:,1:jbeg),2),1,size(data,2));
    else
        for chani=1:size(data,1)
            data(chani,:)=data(chani,:)-median(data(chani,1:jbeg),2);
        end
    end
    
    diary('HBlog.txt')
    warning(['jump? what''','s the noise at ',num2str(jbeg./sRate),'s? zeroed from ',num2str(bads(1)/sRate),' to ',num2str(bads(end)/sRate)]);
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
    meanMEG=double(mean(data));
    meanMEGf = myFilt(meanMEG,BandPassFilt);
    meanMEGf=meanMEGf-median(meanMEGf);
end

%% peak detection on MCG (or ECG) signal
[peaks, Ipeaks]=findPeaks(meanMEGf,1.5,round(sRate*minPeriod)); % 450ms interval minimum
% test if, by chance, the HB field is mainly negative
posHB=true;
if isempty(ECG)
        [peaksNeg, IpeaksNeg]=findPeaks(-meanMEGf,1.5,round(sRate*minPeriod));
    if median(peaksNeg)/median(peaks)>1.1
        diary('HBlog.txt')
        warning('NEGATIVE HB FIELD? if not, average the MEG and give it as ECG');
        diary off
        Ipeaks=IpeaksNeg;
        peaks=-peaksNeg;
        posHB=false;
        %meanMEGf=-meanMEGf;
    end
end
time=1/sRate:1/sRate:size(data,2)/sRate;
if figs
    figure;
    plot(time,meanMEGf)
    hold on
    plot(time(Ipeaks), peaks,'ro')
    title('peak detection on mean MEG (or ECG) trace, OK if many of them are not HB')
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
        if strcmp(cfg.layout,'neuromag306mag.lay')
            [~,magi]=ismember({'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641'},topo.label);
            %topo.avg=topo.avg(chi);
            %topo.label=topo.label(chi);
            cfg.xlim=[1,1];
            cfg.zlim=[-max(abs(topo.avg(magi))) max(abs(topo.avg(magi)))];
            figure;
            ft_topoplotER(cfg,topo);
            title('MAGNETOMETERS, TOPOGRAPHY OF R')
            % cfg.layout='neuromag306planar.lay';
            % grd=topo.avg;
            % grd(chi)=0;
            % cfg.zlim=[-max(abs(grd)) max(abs(grd))];
            % figure;
            % ft_topoplotER(cfg,topo);
            % title('GRADIOMETERS')
        else
            %cfg.channel={'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641'};
            %cfg.interpolation='linear';
            cfg.xlim=[1,1];
            cfg.zlim=[-max(abs(topo.avg)) max(abs(topo.avg))];
            figure;
            ft_topoplotER(cfg,topo);
            title('TOPOGRAPHY OF R')
        end
    else
        warning('no topoplot without labels and layout fields! see figOptions options')
    end
end
topoTrace=median(data(:,Ipeaks),2)'*data;
topoTrace=myFilt(topoTrace,BandPassFilt);
topoTrace=topoTrace-median(topoTrace);
if ~posHB
    meanMEGN=meanMEGf./max(-meanMEGf(1:round(sRate*10)));
    topoTraceN=-topoTrace./max(topoTrace(1:round(sRate*10)));
else
    topoTraceN=topoTrace./max(topoTrace(1:round(sRate*10)));
    meanMEGN=meanMEGf./max(meanMEGf(1:round(sRate*10)));
end




%% check if topo of every peak is correlated to average topo
r=corr(data(:,Ipeaks),median(data(:,Ipeaks),2));
if figs
    figure;
    plot(time,topoTraceN)
    hold on
    plot(time,meanMEGN,'r')
    plot(time(Ipeaks(r>0.5)),meanMEGN(Ipeaks(r>0.5)),'g.');
    legend('topoTrace','meanMEG','r data-topo > 0.5')
end
%% average good HB and make a template
IpeaksR=Ipeaks(r>0.5);
IpeaksR=IpeaksR(IpeaksR>sampBefore);
IpeaksR=IpeaksR(IpeaksR<(size(data,2)-sampBefore));
perEstimate=diff(IpeaksR)/sRate; % estimate period
perEstimate=median(perEstimate(perEstimate<2)); % less than 2s


if posHB
    [temp1e,period1]=makeTempHB(meanMEGf,sRate,IpeaksR,perEstimate,sampBefore,figs);
else
    [temp1e,period1]=makeTempHB(-meanMEGf,sRate,IpeaksR,perEstimate,sampBefore,figs);
end
%% find xcorr between template and meanMEG
if posHB
    xcrPad=XCORR(meanMEGf,temp1e);
else
    xcrPad=XCORR(-meanMEGf,temp1e);
end
% xcr=xcr(lags>=0);
%[~,tempMax]=max(temp1e);
% xcrPad=zeros(size(meanMEGf));
% xcrPad(tempMax:end)=xcr(1:end-tempMax+1);
if figs
    figure;
    plot(time,topoTraceN)
    hold on
    plot(time,meanMEGN,'r')
    if posHB
        plot(time,xcrPad/max(xcrPad),'g');
    else
        plot(time,-xcrPad/max(-xcrPad),'g');
    end
    legend('topoTrace','meanMEG','temp xcorr')
end
%% second sweep
%% find peaks on xcorr trace
[peaks2, Ipeaks2]=findPeaks(xcrPad,1.5,round(sRate*period1*0.65)); % no peaks closer than 60% of period
if figs
    figure;
    if posHB
        plot(time,xcrPad)
        hold on
        plot(time(Ipeaks2), peaks2,'ro')
    else
        plot(time,-xcrPad)
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
% make mcg trace for meanMEG
if posHB
    [temp2e,period2]=makeTempHB(meanMEG,sRate,Ipeaks2in,period1,sampBefore,figs);
else
    [temp2e,period2]=makeTempHB(-meanMEG,sRate,Ipeaks2in,period1,sampBefore,figs);
end

%% test R amplitude
% meanMEGdt=detrend(meanMEG,'linear',round(sRate:sRate:length(meanMEG)));
HighPassSpecObj=fdesign.highpass('Fst,Fp,Ast,Ap',1,2,60,1,sRate);%
HighPassFilt=design(HighPassSpecObj ,'butter');
meanMEGhp = myFilt(meanMEG,HighPassFilt);
% baseline correction again, just in case
meanMEGhp=meanMEGhp-median(meanMEGhp);

[~,maxi]=max(temp2e(1:round(length(temp2e/2))));
bef=find(fliplr(temp2e(1:maxi))<0,1)-1;
aft=find(temp2e(maxi:end)<0,1)-1;
Rlims=[maxi-bef,maxi+aft]; % check where R pulls the template above zero
for HBi=1:length(Ipeaks2in);
    s0=Ipeaks2in(HBi)-bef;
    s1=Ipeaks2in(HBi)+aft;
    x=temp2e(Rlims(1):Rlims(2));
    y=meanMEGhp(s0:s1);
    scalef=-round(log10(max([x,y]))); % scaling factor for polyfit not to complain
    x=x*10^scalef;
    y=y*10^scalef;
    pScaled=polyfit(x,y,1);
    pScaled(2)=pScaled(2).*10^-scalef;
    p(HBi,1:2)=pScaled; %#ok<AGROW>
end
if posHB
    negp=find(p(:,1)<0);
else
    negp=find(p(:,1)>0);
end
if ~isempty(negp)
    p(negp,1:2)=0;
    diary('HBlog.txt')
    warning(['did not get good fit for amplitude test, assume average HB amplitude at ',num2str(Ipeaks2in(negp)/sRate)])
    diary off
end
ampMMfit=p(:,1)+(p(:,2)./temp2e(maxi));
if posHB
    MCG=makeMCG(temp2e,Rlims,Ipeaks2in,ampMMfit,length(meanMEG));
else
    MCG=makeMCG(temp2e,Rlims,Ipeaks2in,-ampMMfit,length(meanMEG));
    MCG=-MCG;
end

%% remove mcg from each chan
%make temp per chan
sBef=maxi-1;
sAft=length(temp2e)-maxi;
HBtemp=HBbyChan(data,sRate,Ipeaks2in,sBef,sAft);
% check SNR per chan
sm50=round(sRate*0.05);
s0=maxi-sm50*3;
s1=maxi-sm50;
s2=maxi+sm50;
n=std(HBtemp(:,s0:s1)'); %#ok<*UDIM>
s=std(HBtemp(:,s1:s2)');
snr=s./n;
badSNR=find(snr<=snrThr);
if figs
    timeTemp=1/sRate:1/sRate:size(HBtemp,2)/sRate;
    timeTemp=timeTemp-maxi/sRate;
    figure;plot(timeTemp,HBtemp','k');hold on;
    %plot(maxi,HBtemp(find(HBsnr>1.1),maxi),'g.')
    [~,minI]=min(snr);
    plot(timeTemp,HBtemp(minI,:),'r');
    if isempty(badSNR) || snrThr==0
        title('HB for all channels,red trace is worst SNR channel')
    else
        plot(timeTemp(maxi),HBtemp(badSNR,maxi),'r.')
        title(['HB for all channels. red trace is worst SNR channel, dots mark channels with SNR < ',num2str(snrThr)])
    end
end
if ~isempty(badSNR)
    if length(find(snr<=snrThr))>size(data,1)/2
        error('too many channels have poor SNR. was there artifact?')
    end
    if ~snrThr==0
        HBtemp(badSNR,:)=0;
        diary('HBlog.txt')
        disp(['not including bad SNR channels, index no. ',num2str(badSNR)])
        diary off
    end
end
if posHB
    MCGall=makeMCGbyCh(HBtemp,Rlims,Ipeaks2in,ampMMfit,length(meanMEG),maxi);
else
    MCGall=makeMCGbyCh(HBtemp,Rlims,Ipeaks2in,-ampMMfit,length(meanMEG),maxi);
end
cleanData=data-MCGall;
figure;
plot(time,MCG,'k')
hold on
plot(time,mean(data),'r')
plot(time,mean(cleanData),'g')
legend('MCG from template', 'mean MEG','mean clean MEG')

Rtopo=HBtemp(:,maxi);
if figs
    topo.avg=Rtopo;
    cfg.xlim=[1,1];
    if strcmp(cfg.layout,'neuromag306mag.lay')
        cfg.zlim=[-max(abs(topo.avg(magi))) max(abs(topo.avg(magi)))];
        figure;
        ft_topoplotER(cfg,topo);
        title('MAGNETOMETERS, TOPOGRAPHY OF R, 2nd sweep')
    else
        cfg.zlim=[-max(abs(Rtopo)) max(abs(Rtopo))];
        figure;
        ft_topoplotER(cfg,topo);
        title ('TOPOGRAPHY OF R, 2nd sweep')
    end
end
display(['HB period (2nd sweep) is ',num2str(period2),'s']);
function [tempe,period2]=makeTempHB(trace,sRate,peakIndex,period,sampBefore,figs)
betweenHBs=0.7; % after the T wave, before next qrs, 0.7 of period
HB=zeros(size(trace,1),sampBefore*2+1);
for epochi=1:length(peakIndex)
    HB=HB+trace(:,peakIndex(epochi)-sampBefore:peakIndex(epochi)+sampBefore);
end
HB=HB/epochi;
period2=[]; %#ok<NASGU>
HBxc=xcorr(HB);
[~,ipxc]=findPeaks(HBxc,1.5,sRate*period*0.6); % index peak xcorr
if length(ipxc)>1
    nextPeak=ceil(length(ipxc)/2+0.5);
    period2=(ipxc(nextPeak)-ipxc(nextPeak-1))/sRate;
    % else
    xcorrByseg=false;
    if length(trace)>=2^20 % test if version is later than 1011b
        ver=version('-release');
        if ~strcmp(ver,'2011b')
            if str2num(ver(1:end-1))<=2011
                xcorrByseg=true;
            end
        end
    end
    if xcorrByseg
        trace1=trace;
        difs=[];
        while length(trace1)>length(HB)*3%2^20
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
            % FIXME xcorr (fftfilt) won't take it for more than 2^20
        end
        period2=median(difs)/sRate;
    else
        HBxc1=xcorr(trace,HB);
        [~,ipxc]=findPeaks(HBxc1,1.5,sRate*period*0.6); % index peak xcorr
        period2=median(diff(ipxc))/sRate;
    end
else
    warning('could not find cross correlation within extended template, guessing period')
    period2=period;
end
temp=HB(sampBefore-round(sRate*(1-betweenHBs)*period2):sampBefore+round(sRate*betweenHBs*period2));
edgeRepressor=ones(size(temp));
ms20=round(sRate/50);
reducVec=0:1/ms20:1;
reducVec=reducVec(1:end-1);
edgeRepressor(1:length(reducVec))=reducVec;
edgeRepressor(end-length(reducVec)+1:end)=fliplr(reducVec);
tempe=temp-median(temp);
tempe=tempe.*edgeRepressor;
time=1/sRate:1/sRate:length(temp)/sRate;
time=time-(1-betweenHBs)*period2;
if figs
    figure;
    plot(time,tempe,'g')
    title('template HB')
end
tempe=double(tempe);

function MCG=makeMCG(temp,Rlims,Ipeaks,amp,lengt)
MCG=zeros(1,lengt);
[~,maxTemp]=max(temp(1:round(length(temp/2))));
for HBi=1:length(Ipeaks);
    s0=Ipeaks(HBi)-maxTemp+1;
    s1=Ipeaks(HBi)+length(temp)-maxTemp;
    MCG(s0:s1)=temp;
    MCG(s0+Rlims(1)-1:s0+Rlims(2)-1)=temp(Rlims(1):Rlims(2))*amp(HBi);
end
function tempe=HBbyChan(trace,sRate,peakIndex,sampBefore,sampAfter)
HB=zeros(size(trace,1),sampBefore+1+sampAfter);
% average HBs
for epochi=1:length(peakIndex)
    HB=HB+trace(:,peakIndex(epochi)-sampBefore:peakIndex(epochi)+sampAfter);
end
HB=HB/epochi;
% reduce edges to vero
edgeRepressor=ones(1,size(HB,2));
ms20=round(sRate/50);
reducVec=0:1/ms20:1;
reducVec=reducVec(1:end-1);
edgeRepressor(1:length(reducVec))=reducVec;
edgeRepressor(end-length(reducVec)+1:end)=fliplr(reducVec);
tempe=HB-repmat(mean(HB(:,[1:ms20,end-ms20:end]),2),1,size(HB,2));
tempe=tempe.*repmat(edgeRepressor,size(HB,1),1);
function MCG=makeMCGbyCh(temp,Rlims,Ipeaks,amp,lengt,maxTemp)
MCG=zeros(size(temp,1),lengt);
for HBi=1:length(Ipeaks);
    s0=Ipeaks(HBi)-maxTemp+1;
    s1=Ipeaks(HBi)+length(temp)-maxTemp;
    MCG(:,s0:s1)=temp;
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
[~,tempMax]=max(y);
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
        xcrPad=zeros(size(trace2));
        xcrPad(tempMax:end)=xcCurrent(1:end-tempMax+1);
        xcr=[xcr,xcrPad];
        % FIXME xcorr (fftfilt) won't take it for more than 2^20
    end
    
else
    [xcr,lags]=xcorr(x,y);
    xcr=xcr(lags>=0);
    xcrPad=zeros(size(x));
    xcrPad(tempMax:end)=xcr(1:end-tempMax+1);
    xcr=xcrPad; % sorry for switching variables
end
% this is for older than 2011b

