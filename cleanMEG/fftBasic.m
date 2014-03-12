function [fourier,freq]=fftBasic(rows,Fs)
% rows = data with rows for channels
% Fs = Sampling frequency

L = size(rows,2)/Fs;                     % Length of signal
NFFT = round(Fs); % this gives bins of roughly  1Hz
for seci=1:NFFT:(size(rows,2)-NFFT)
    if seci==1
        Y = fft(rows',NFFT);
    else
        Y = Y+fft(rows(:,seci:seci+NFFT)',NFFT);
    end
end
Y=Y./length(1:NFFT:(size(rows,2)-NFFT)); % average over seconds
fourier=Y(1:floor(NFFT/2)+1,:);
freq = Fs/2*linspace(0,1,NFFT/2+1);
fourier=fourier';
freq=freq(2:end);
fourier=fourier(:,2:end);
end