function [peaks, Ipeaks] = findPeaks(x,th, deadT)
% find all peaks in x
% [peaks, Ipeaks] = findPeaks(x,th);
%
% x      - column vector
% th     - consider only peaks which are aboveth*std
% deadT  - minimal number of samples between succesive peaks.  When closer
%          choose the bigest. [default 0];
% peaks  - values of peaks
% Ipeaks - indices in x where the peaks are

% nov-2007 MA
%  UPDATES
% Jan-2010 mean and std are ignoring NaNs.  MA
% Jul-2010 Default deadt is 0.  MA

if ~exist('th','var'), th=0; end
if isempty(th), th=0; end
if ~exist('deadT','var'), deadT=0; end
if isempty(deadT), deadT=0; end
[r,c] = size(x);
if (r<c), x=x'; end
dx = diff(x);
d2x = diff(dx)<0;
changeDir = dx(1:end-1).*dx(2:end)<0;
Ipeaks = find(changeDir&d2x)+1;
peaks = x(Ipeaks);
if th>0
    sd =stdWnan(x);
    mn=meanWnan(x);
    smallP = find(peaks<(mn+th*sd));
    peaks(smallP)=[];
    Ipeaks(smallP)=[];
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
