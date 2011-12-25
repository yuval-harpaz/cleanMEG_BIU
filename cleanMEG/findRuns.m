function [runLength, nameId] = findRuns(I)
% find runs of length 1,2,... with the same value in I
%
%  I - vector of values
%  runLength -  list of runs-lengthes as they appeared in I
%  nameId    -  the value in I with the appropriate  run Length
%  e.g  I=[1,2,3,3,1,4,1,1,1,2], will produce
%  runLength    nameId
%     1            1
%     1            2
%     2            3
%     1            1
%     1            4
%     3            1
%     1            2

%  Sep-2005  MA

lI = length(I);

k1=2;      % index into I
kL=1;      % length of run
ii = 1;    % serial number of run
i1 = I(1); % name of current value in I
while(k1<=lI)
    if I(k1)~=i1  % a new run starts
       runLength(ii) = kL;
       nameId(ii)    = i1;
       ii = ii+1;
       kL=1;
       i1 = I(k1);
   else  % in the same 
       kL = kL+1;
   end
   k1 = k1+1;
end
% add the last one
runLength(ii) = kL;
nameId(ii)    = i1;

return