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
cleanPeriods=findCleanPeriods(fileName);
goodPeriods=sumGoodPeriods(fileName,cleanPeriods,[]);

for segi=1:size(goodPeriods,2)
    p=pdf4D(fileName);
    cleanCoefs = createCleanFile(p, fileName,...
        'byLF',256 ,'Method','Adaptive',...
        'xClean',[4,5,6],...
        'CleanPartOnly',[goodPeriods(1,segi) goodPeriods(2,segi)],...
        'chans2ignore',[74,204],...
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
