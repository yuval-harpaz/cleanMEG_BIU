function [fourier,freq]=fftRaw(fileName,startTime,duration,keepSegments)
% fileName - 4D pdf file 'c,rf...'
% startTime and duraton can be given in sec, not that parts of seconds will
% not be ffted in the end of the data.
% keep segments for sec by sec results to be kept unaveraged
% similar to fftBasic but I added baseline correction for DC data and no
% short (less than 1sec) data. 

if ~exist('fileName','var')
    try
        fileName=source;
    end
end
if ~exist('startTime','var')
    startTime=0;
end
if ~exist('keepSegments','var')
    keepSegments=false;
end
p=pdf4D(fileName);
chi = channel_index(p, 'meg', 'name');
Fs=get(p,'dr');
startSamp=max(round(startTime.*Fs), 1);
if ~exist('duration','var')
    hdr=get(p,'header');
    endTime=hdr.epoch_data{1}.epoch_duration;
    endSamp=ceil(endTime.*Fs);
else
    endSamp=ceil(duration.*Fs)+startSamp;
end
display(['reading ',fileName]);
rows = read_data_block(p,[startSamp endSamp],chi);
    

L = size(rows,2)/Fs;                     % Length of signal
NFFT = round(Fs); % this gives bins of roughly  1Hz
secCount=0;
if NFFT<size(rows,2)
    for seci=1:NFFT:(size(rows,2))
        try % should fail for too short, end of the rows
            dataSec=rows(:,seci:seci+NFFT)';
            secCount=secCount+1;
            dataSec=dataSec-repmat(mean(dataSec),size(dataSec,1),1);
            if seci==1
                Y = fft(dataSec,NFFT);
            else
                if keepSegments
                    Y(:,:,secCount)=fft(dataSec,NFFT);
                else
                    Y = Y+fft(dataSec,NFFT);
                end
            end
            
        end
    end
    if ~keepSegments
        Y=Y./secCount; % average over seconds
    end
else % short data, less than a second
    error('short data?')
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