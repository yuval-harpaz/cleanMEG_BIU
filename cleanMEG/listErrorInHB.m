function tList = listErrorInHB(cleanCoefs)
% list the times at which there was a difficulty with finding Hear Beats
%    tList = listErrorInHB(cleanCoefs);
% cleanCoefs - array of structures with info on heart beats, as generated
%              by createCleanFileE
% tList      - list of times in seconds where some difficulty was
%              found

% Jan-2010  MA

%% initialize
numPieces = length(cleanCoefs);
lastTime = cleanCoefs(numPieces).HBparams.lastT;
numProblems = round(lastTime);
tList = zeros(1, numProblems);
pIndex =1;

%% find errors
for pNo=1:numPieces
    err = cleanCoefs(pNo).HBparams.Errors;
    T0 = cleanCoefs(pNo).HBparams.firstT;
    if isfield(err,'shortHB')
        t = err.shortHB;
        for ii = 1:length(t)
            tList(pIndex) = T0 +t(ii);
            pIndex = pIndex+1;
        end
    end
    if isfield(err,'longHB')
        t = err.longHB;
        for ii = 1:length(t)
            tList(pIndex) = T0 +t(ii);
            pIndex = pIndex+1;
        end
    end
    if isfield(err,'smallHB')
        t = err.smallHB;
        for ii = 1:length(t)
            tList(pIndex) = T0 +t(ii);
            pIndex = pIndex+1;
        end
    end
    if isfield(err,'bigHB')
        t = err.bigHB;
        for ii = 1:length(t)
            tList(pIndex) = T0 +t(ii);
            pIndex = pIndex+1;
        end
    end
end  % end of going over all pieces

%% wrapUp
tList(pIndex:end)=[];
tList = sort(tList);

return
