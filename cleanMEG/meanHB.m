function [avgHB,times]=meanHB(data,sRate,HBtimes,plt);
if ~exist ('plt','var')
    plt=false;
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
    