function whichLF = findLFbit(trig,samplingRate)
% find which bit oscillates in ~ line frequency

% Jan 2009  MA

%% initialize
lf1 = samplingRate/50;
lf2 = samplingRate/60;
bit = uint16(1);
trig= uint16(trig);

%% search for the bit
for ii = 1:16
    testTrig = bitand(trig, bit);
    whereUp=find(diff(double(testTrig==bit))>0);
    if ~isempty(whereUp)
        meanLF = mean(diff(whereUp));
        if (lf1-1)<meanLF && meanLF<(lf1+1)
            whichLF = single(bit);
            return
        elseif (lf2-1)<meanLF && meanLF<(lf2+1)
            whichLF = single(bit);
            return
        end
    end
    bit = 2*bit;
end

whichLF = [];
return