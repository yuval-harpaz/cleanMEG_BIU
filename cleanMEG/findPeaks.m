function [peaks, Ipeaks] = findPeaks(x,th, deadT)
% find all peaks in x
% [peaks, Ipeaks] = findPeaks(x,th,deadT);
%
% x      - column vector
% th     - consider only peaks which are above  th*std. when th is text it can also be
% percentile (e.g. '%95') or a fixed value (e.g. 'v0.0007')
% deadT  - minimal number of samples between succesive peaks.  When closer
%          choose the bigest. [default 0];
% peaks  - values of peaks
% Ipeaks - indices in x where the peaks are

% nov-2007 MA
%  UPDATES
% Jan-2010 mean and std are ignoring NaNs.  MA
% Jul-2010 Default deadt is 0.  MA
% Mar-2011 Peak of plataus are also peaks!!(were ignored previousl!) MA

%% initialize
if ~exist('th','var'), th=0; end
if isempty(th), th=0; end
if ~exist('deadT','var'), deadT=0; end
if isempty(deadT), deadT=0; end
[r,c] = size(x);
if (r<c), x=x'; end

%% find peaks and plataus
dx = diff(x);
% d2x = diff(dx)<0;
% changeDir = dx(1:end-1).*dx(2:end)<0;
% Ipeaks = find(changeDir&d2x)+1;

% new approach to detect also plataus
J=find(dx~=0)';
P=dx(J)>0;
N=dx(J)<0;
% K = find(P(1:end-1) & N(2:end));
Ipeaks = J(P(1:end-1) & N(2:end)) +1;
% if at the peaks there are several identical values choose the middle
for ii = 1:length(Ipeaks)
    platau = 0;
    thisP = x(Ipeaks(ii));
    thisI = Ipeaks(ii);
    for jj = thisI:length(Ipeaks)
        if x(jj)==thisP
            platau = platau+1;
        else
            break
        end
    end
    Ipeaks(ii) = thisI +floor(platau/2);
end

%% lint the peaks according to size and dead time
peaks = x(Ipeaks);
if ischar(th)
    if strcmp(th(1),'%') % percentile
        percentile = prctile(x,str2double(th(2:end)));
        smallP = find(peaks<percentile);
        peaks(smallP)=[];
        Ipeaks(smallP)=[];
    elseif strcmp(th(1),'v') % fixed value, in volt or whatever unit x has
        smallP = find(peaks<str2num(th(2:end)));
        peaks(smallP)=[];
        Ipeaks(smallP)=[];
    else
        error ('text threshold can only start with % or v')
    end
elseif th>0
    sd =stdWnan(x);
    mn=meanWnan(x);
    smallP = find(peaks<(mn+th*sd));
    peaks(smallP)=[];
    Ipeaks(smallP)=[];
end
if isempty(Ipeaks)
    return
end
if deadT>1
    eliminate = true;  % initial start
    while sum(eliminate>0)
        eliminate = false(size(Ipeaks));
        ii = 2;
        oldI = Ipeaks(1);
        oldP = peaks(1);
        while ii<=length(Ipeaks)
            if (Ipeaks(ii)-oldI)<=deadT
                if peaks(ii)>oldP
                    eliminate(ii-1)=true;
                    oldI = Ipeaks(ii);
                    oldP = peaks(ii);
                else
                    eliminate(ii)=true;
                end
            else
                oldI = Ipeaks(ii);
                oldP = peaks(ii);
            end
            ii = ii+1;
        end
        peaks(eliminate)=[];
        Ipeaks(eliminate)=[];
    end
end
return
