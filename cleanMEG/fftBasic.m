function [fourier,freq]=fftBasic(rows,Fs)
% rows = data with rows for channels
% Fs = Sampling frequency
L = size(rows,2)/Fs;                     % Length of signal
NFFT = round(Fs); % this gives bins of roughly  1Hz
Y = fft(rows',NFFT);
fourier=Y(1:floor(NFFT/2)+1,:);
freq = Fs/2*linspace(0,1,NFFT/2+1);
fourier=fourier';
freq=freq(2:end);
fourier=fourier(:,2:end);
end