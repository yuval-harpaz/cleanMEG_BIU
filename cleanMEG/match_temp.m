function [snr,signal]=match_temp(dataVec,temp,tempZero)
% here you look for a match between a template and data. tempZero is the
% time of interest in the template (in samples), like the maximum point of the event.

numSlides=length(dataVec)-length(temp)+1; % how many times the template will be shifted to match the data
if numSlides<1
    error('the data has to be larger than the template')
end
% baseline correction (sum = 0)
tempBlc = temp-mean(temp);
% normalization (sum of squares = 1)
tempBlcNorm=tempBlc./sqrt(sum(tempBlc.*tempBlc));
% try
%     matlabpool;
% end
snr=zeros(size(dataVec));
signal=snr;
for slidei=1:numSlides
    data=dataVec(slidei:slidei+length(temp)-1);
    data=data-mean(data);
    data=data./sqrt(sum(data.*data));
    Proj = sum(data.*tempBlcNorm);
    Signal = Proj^2;
    %Total = sum(data.*data);
    Total = 1;
    Error = Total - Signal;
    if Error<1e-14
        warning(['template matches the signal too well at lag ',num2str(slidei)])
        if Error==0
            Error=1e-14;
            display('replacing Error - zero with 1e-14')
        end
    end
    SNR=Signal/Error;
    Sign=Proj/abs(Proj);
    signal(tempZero+slidei-1)=Sign*Signal;
    snr(tempZero+slidei-1)=Sign*SNR;
end



