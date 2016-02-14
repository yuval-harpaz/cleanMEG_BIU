function [kursum,kurmax,rep,kurByWin]=kurLoop(X,Nsamp,overlap,flag)
% computes kurtosis with a sliding window
% X and flag are identical to 'kurtosis' function, but note the dimentions:
% a data matrix (or a vector), raws for channels, columns for time samples.
% Nsamp is how many samples to include in a window
% overlap should  0.5 or 0
if ~exist('overlap','var')
    overlap=0;
end
if ~exist('flag','var')
    flag=[];
end

kursum=zeros(size(X,1),1);
kurmax=kursum;
samp0=round(1:(Nsamp-overlap*Nsamp):size(X,2));
samp1=samp0+round(Nsamp)-1;
out=find(samp1>size(X,2),1);
if ~ isempty(out)
    samp0=samp0(1:out-1);
    samp1=samp1(1:out-1);
end
rep=length(samp0);
kurByWin=[];
for wini=1:rep
    kurwin=kurtosis(X(:,samp0(wini):samp1(wini)),flag,2);
    % kurwin(kurwin<0)=0;
    kursum=kursum+kurwin;
    kurmax=max([kurmax,kurwin],[],2);
    if nargout==4
        kurByWin=[kurByWin,kurwin];
    end   
end


