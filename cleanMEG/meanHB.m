function [avgHB,times]=meanHB(data,sRate,HBtimes,plt);
if ~exist ('plt','var')
    plt=false;
end
if ischar(data)
    if isempty(which('ft_preprocessing'))
        error('data file name needs FieldTrip in path')
    else
        cfg=[];
        cfg.dataset=data;
        cfg.demean='yes';
        cfg.continuous='yes';
        cfg.channel='MEG';
        data=ft_preprocessing(cfg);
        if isempty(sRate)
            sRate=data.fsample;
        end
        data=data.trial{1};
    end
end
HBsamp=round(HBtimes*sRate);
avgHB=zeros(size(data,1),round(sRate)*2+1);
counter=0;
for HBi=1:length(HBsamp)
    try % it will fail in the end
        avgHB=avgHB+data(:,HBsamp(HBi)-round(sRate):HBsamp(HBi)+round(sRate));
        counter=counter+1;
    end
end
times=-round(sRate):round(sRate);
times=times/sRate;
avgHB=avgHB./counter;
avgHB=avgHB-repmat(median(avgHB,2),1,size(avgHB,2));
if plt
    figure;plot(times,avgHB,'k');
end
    