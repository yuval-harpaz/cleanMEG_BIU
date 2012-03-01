function [goodPeriods,goodsamples1]=sumGoodPeriods(fileName,cleanPeriods,chanNumThr)
% creates a list of good periods which is good for all channels
% cleanPeriods is created by findCleanPeriods
% fileName='c,rfhp0.1Hz';
% chunNumThr is the threshold of bad channels above which a segment is
% taged as bad (defalt=20).
if ~exist('chanNumThr','var')
    chanNumThr=[];
end
if isempty(chanNumThr)
    chanNumThr=30;
end
p=pdf4D(fileName);
hdr=get(p,'header');
sampleRate=get(p,'dr');
goodsamples1=zeros(1,hdr.epoch_data{1,1}.pts_in_epoch);
for chani=1:length(cleanPeriods)
    for periodi=1:size(cleanPeriods{1,chani},2)
        endtrunk=0;
        if periodi==size(cleanPeriods{1,chani},2)
            endtrunk=-1;
        end
        goodsamples1(ceil(cleanPeriods{1,chani}(1,periodi)*sampleRate):...
            floor(cleanPeriods{1,chani}(2,periodi)*sampleRate)+endtrunk)=...
            goodsamples1(ceil(cleanPeriods{1,chani}(1,periodi)*sampleRate):...
            floor(cleanPeriods{1,chani}(2,periodi)*sampleRate)+endtrunk)...
            +1;
    end
end
goodsamp=goodsamples1>=(length(cleanPeriods)-chanNumThr);
figure;plot(goodsamples1,'b');hold on;plot((length(cleanPeriods)-chanNumThr)*goodsamp,'k');
legend('NUMBER OF GOOD CHANNELS PER SAMPLE','GOOD EPOCHS AFTER APPLYING THRESHOLD')
diffgood=diff(goodsamp);
goodon=find(diffgood==1);
goodoff=find(diffgood==-1);
if goodon(1)>goodoff(1) % beggining of file is good
    goodon(1,2:end+1)=goodon;
    goodon(1,1)=1;
end
if goodon(end)>goodoff(end) % end of file is good
    goodoff(end+1)=hdr.epoch_data{1,1}.pts_in_epoch;
end
if ~length(goodon)==length(goodoff)
    error('size mismatch, sort it out, sorry');
end
goodPeriods=goodon/sampleRate;
goodPeriods(2,:)=goodoff/sampleRate;
end


    