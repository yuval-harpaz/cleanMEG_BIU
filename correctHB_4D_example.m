% Example of cleaning 4D data of HeartBeat artifact
% open newData and rewrite_pdf for more options if the default doesn't work
% better clean LF before HB

% clean HB
newData=correctHB;
% write it to a new file
rewrite_pdf(newData);
% rename the rw_* file if you want to


%% BIU users
% better clean line frequency and accelerometers first
% use createCleanFile with 'HeartBeat',0, like this
fileName = 'c,rfhp0.1Hz';
p=pdf4D(fileName);
cleanCoefs = createCleanFile(p, fileName,'byLF',256 ,'Method','Adaptive','xClean',[4,5,6],'byFFT',0,'HeartBeat',0);
% then you get a file called 'xc,lf_c,rfhp0.1Hz' or something similar
% then clean the heartbeat like this (it looks for xc,lf_* file)
newData=correctHB; % or newData=correctHB([],[],1); if you want lots of plots
rewrite_pdf(newData); % or rewrite_pdf(newData,[],[],'xc,lf,hb') if you want standard BIU naming



