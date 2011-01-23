function  bin = gaussBin(sigma,extent,hollow)
% bin = gaussBin(sigma,extent,hollow) to create a Guassian bin with 
% Sd of sigma having 2*extent+1 bins and with center bin multiplied by
% hollow.  The whole final bin is normalized to area of 1

%  MA
% UPDATES
% Aug-2008  HOLLOW is added
%% initialize
if ~exist('hollow', 'var'), hollow=[]; end
if isempty(hollow), hollow=1; end

if sigma == 0;
   bin = 1;
   return
end

if ~exist('extent', 'var')
   sz = ceil(3.5*sigma);
else
   sz=ceil(extent*sigma);
end
var = sigma*sigma;

%% compute gaussian
strt = round(0-sz);
ends = round(0+sz);
ng = strt:1:ends;
mid = find(ng==0);
M = (ng.*ng)/(2*var);
G = exp(-M);

%% make a hollow bin and normalize
G(mid) = hollow*G(mid);
bin = G/sum(G);

return
