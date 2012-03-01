% Avoiding muscle artifact periods, cleaning a file segment by segment
% it is possible to clean good periods by 'CleanPartOnly'. however, after
% the cleaning the beggining of the file tryClean function doesn't work.
% for this purpose the first run of tryClean saves an output file
% 'tryCleanOP.mat' which is used for the consequtive periods.
% for automatic detection of muscle artifact use findCleanPeriods. this
% gives beggining and end times of good periods for every channel. to
% decide which periods are really good for most of the channels use
% sumGoodPeriods. the 3rd input arg is a threshold- the number of bad
% channels per time point to define this timepoint as bad (default is 30).
%
% here is an example, copy and put it elsewhere before modifying it.


fileName='c,rfhp0.1Hz';
fileNameOrig=fileName;

% finding good periods for every channels. based on fft run with a 
% window of 2s sliding in steps of 0.5s. gives for every channel collumns
% of beginning and end times of clean periods.
cleanPeriodsAllChans=findCleanPeriods(fileName);

% deciding which time points are realy clean. strictest is when for a given
% time point there is no bad channel. for this the third argument
% (chanNumThr) has to be 1. the default is 20, that is if 20 channels or more are
% noisy at a sertain a time point it is considered as bad.
cleanPeriods=sumGoodPeriods(fileName,cleanPeriodsAllChans,[]);
% removing too short good periods, less than 5s long
notTooShort=find((cleanPeriods(2,:)-cleanPeriods(1,:))>=5);
cleanPeriods=cleanPeriods(:,notTooShort);
save cleanPeriods cleanPeriods
for segi=1:size(cleanPeriods,2)
    p=pdf4D(fileName);
    cleanCoefs = createCleanFile(p, fileName,...
        'byLF',512 ,'Method','Adaptive',...
        'xClean',[4,5,6],...
        'CleanPartOnly',[cleanPeriods(1,segi) cleanPeriods(2,segi)],...
        'outFile','temp2',...
        'noQuestions',1,...
        'byFFT',0,...
        'HeartBeat',[],... % for automatic HB cleaning change 0 to []
        'maskTrigBits', 512);
    if exist('temp1','file')
        !rm temp1
    end
    !mv temp2 temp1
    fileName='temp1';
    %close all
end
eval(['!mv temp1 per_',fileNameOrig]);
p=pdf4D(fileNameOrig);
hdr=get(p,'header');
lat=[1 hdr.epoch_data{1,1}.pts_in_epoch];
chi=channel_index(p,'MEG','name');
orig=mean(read_data_block(p,lat,chi),1);
p=pdf4D(['per_',fileNameOrig]);
clean=mean(read_data_block(p,lat,chi),1);
figure;plot(orig,'r');hold on;plot(clean,'g');
title('AVERAGED CHANNELS')
legend ('OLD','NEW')
