function [whichLF, OKbit] = findLFbit(trig,samplingRate)
% find which bits oscillates in ~ line frequency

% Jan 2009  MA
%      UPDATES
% Jul-2011  All bits which are near 50 and 60Hz are returned  MA

%% initialize
lf1 = samplingRate/50;
lf2 = samplingRate/60;
bit = uint16(1);
trig= uint16(trig);
OKbit = false(1,16); 

%% search for the bit
for ii = 1:16
    testTrig = bitand(trig, bit);
    whereUp=find(diff(double(testTrig==bit))>0);
    if ~isempty(whereUp)
        meanLF = mean(diff(whereUp));
        if (lf1-1)<meanLF && meanLF<(lf1+1)
            OKbit(ii) = true;
%             whichLF = single(bit);
%             return
        elseif (lf2-1)<meanLF && meanLF<(lf2+1)
            OKbit(ii) = true;
%             whichLF = single(bit);
%             return
        end
    end
    bit = 2*bit;
end

%% wrap up
if sum(OKbit) == 0  % none found
    whichLF = [];
else
    firstOK = find(OKbit,1)-1;
    whichLF = uint16(round(2^firstOK));
end
return