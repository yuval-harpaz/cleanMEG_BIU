function [chnSorted, chiSorted] = sortMEGnames(chn,chi)

%Aug-2008  MA
if size(chn,2)==1
    chn=chn';
end
numChannels = size(chn,2);
chnSorted=cell(1,numChannels);
for ii=1:numChannels
    name = chn{ii};
    num2pad = 4-length(name);
    switch num2pad
        case 1
            name = [name(1) '0' name(2:end)];
        case 2
            name = [name(1) '00' name(2:end)];
    end
    chnSorted{ii}=name;
end

[chnSorted, I] = sort(chnSorted);
chiSorted = chi(I);  %I; %(I);

return