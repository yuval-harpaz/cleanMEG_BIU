function outTrig=clearBits(inTrig, bits2clear)
% clear the bits of inTrig which are in bits2clear
% works on 16 bits only

% Nov-2008 MA

% find what kind of variable is inTrig
A=whos('inTrig');
kind = A.class;
%mask the bits
inTrig = uint16(inTrig);
bits2clear = bitcmp(uint16(bits2clear));
outTrig = bitand(inTrig, bits2clear);
% convert to the input type
eval(['outTrig = ' kind '(outTrig);'])
return
