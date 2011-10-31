function s = stdWnan(x,dim)
% compute the standard deviation even if there are NaNs along dimension d
% s = stdWnan(x,dim);
%
% x   - the data, either array up to 4D matrix
% dim - dimension along which to compute, if
%       x  is a singleton the std along its greater then 1 dimension is used.
% s   - the standard deviation normalized by N-1, where N is the number of non NaN
%       elements.  If N is 1 or 0 a NaN is returned

% Jan 2008 MA
% UPDATES
%  aug-2009 extende for 4D

%% initialize
if ~exist('dim','var'), dim=1; end
if isempty(dim), dim=1; end
%% compute
s=sqrt(varWnan(x,dim));
return

