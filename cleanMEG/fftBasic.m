function [fourier,freq]=fftBasic(rows,Fs,keepSegments)
% rows = data with rows for channels
% Fs = Sampling frequency
% you can give 4D filename too as first argument.
if ischar(rows)
    fileName=rows;
    p=pdf4D(fileName);
    chi = channel_index(p, 'meg', 'name');
    display(['reading ',fileName]);
    rows = read_data_block(p,[],chi);
    Fs=get(p,'dr');
end
if ~exist('keepSegments','var')
    keepSegments=false;
end
L = size(rows,2)/Fs;                     % Length of signal
NFFT = round(Fs); % this gives bins of roughly  1Hz
secCount=1;
if NFFT<size(rows,2)
for seci=1:NFFT:(size(rows,2))
    if seci==1
        Y = fft(rows',NFFT);
    else
        try % should fail for too short, end of the rows
        secCount=secCount+1;
        if keepSegments
            Y(:,:,secCount)=fft(rows(:,seci:seci+NFFT)',NFFT);
        else
            Y = Y+fft(rows(:,seci:seci+NFFT)',NFFT);
        end
        end
    end
end
if ~keepSegments
    Y=Y./secCount; % average over seconds
end
else % short data, less than a second
    NFFT=size(rows,2);
    Y=fft(rows',NFFT);
end
fourier=Y(1:floor(NFFT/2)+1,:,:);
freq = Fs/2*linspace(0,1,NFFT/2+1);
if keepSegments
    fourier=permute(fourier,[2,1,3]);
else
    fourier=fourier';
end
freq=freq(2:end);
fourier=fourier(:,2:end,:);
end