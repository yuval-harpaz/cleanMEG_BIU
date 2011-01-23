function [result, varFactor] = myConv(dataA, smoothA, Window, normalize)
% [result, varFactor] = myConv(data, smooth, Window, normalize);
% smooth data with a smoothing bin, expand data symetricallty around edges
% to avoid edge effects, return also a factor for variance reduction.
%
% data      - data to be smoothed
% smooth    - the bin for smoothing
% Window    - multiple every piece by Window before convolving
% normalize - 'NORMALIZE' if to normalize the area under smooth to 1
% result    - the smoothed data
% varFactor - how much the data variance is reduced by smoothing
%
% Assumptions:
%  data and smooth have an odd length
%  adjacant data points are iid random variables
% NOTE smooth is normalizes such that sum(smooth) =1

% see copybook 324, part I, p16
% April 2002, MA
% UPDATES
%  Apr-2005  the longer array is assumed to be the data
%   extending the data is done only by as much as is needed
%  Feb-2008  result is the same length as data 
%  Jan-2010  Window added  MA

%% initialize input params
if ~exist('normalize','var'), normalize='normalize'; end
toNormalize = strcmpi(normalize,'normalize');
if ~exist('Window','var'), Window=[]; end
if isempty(Window)
    toWindow = false;
else
    toWindow = true;
end

%% find which is longer
if length(dataA)>=length(smoothA)
    data = dataA;
    smooth = smoothA;
else
    disp('2-nd argument is longer - assuming it is the data');
    smooth = dataA;
    data = smoothA;
end
[a,b] = size(data);
if a>b
    X=data';
else
    X=data;
end
% make sure smooth is a row vector
if size(smooth,1)>1, smooth = smooth'; end
if toWindow  % make sure length of smoothing kernel and Window are identical
    if size(Window,1)>1, Window = Window'; end
    if size(Window,2)~=size(smooth,2)
        error('MATALB:statistic:ImproperParameters',...
              ' Window and smoothing window should be of equal length');
    end
end
%% convolve
if toNormalize
    normF = sum(smooth);
else
    normF=1;
end
varFactor = sum(smooth.*smooth)/(sum(smooth)*sum(smooth));
middleX = floor(length(X)/2)+1;
% middleS = floor(length(smooth)/2)+1;
ls = length(smooth);
lsHalf = round(ls/2);
ld = length(data);
expandX = [fliplr(X(2:ls+1)),X,fliplr(X(end-ls:end-1))];

if ~toWindow
    convolved = conv(expandX, smooth);
    middleCon = floor(length(convolved)/2) + 1;
    result = convolved(middleCon-middleX+1 : middleCon+middleX-1)/normF;
else % do in parts of Window length
    ls = length(smooth)-1;
    result = nan(1,length(X));
    for ii = 1:length(X)
        iStrt = ii+lsHalf;
        iEnd = iStrt+ls;
        x=expandX(iStrt:iEnd);
        x = x.*Window;
        x = sum(x.*smooth);
        result(ii)=x;
    end
    return
end

%% make sure the same length as DATA
dn = ld-length(result);
if dn>0 % add one point at the end
    result = convolved(middleCon-middleX+1 : middleCon+middleX)/normF;
elseif dn<0 % remove one point
    result = convolved(middleCon-middleX+1 : middleCon+middleX-2)/normF;
end

return
