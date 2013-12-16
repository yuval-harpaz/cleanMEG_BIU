% Example of cleaning 4D data of HeartBeat artifact
% open newData and rewrite_pdf for more options if the default doesn't work


% clean HB
newData=correctHB;
% write it to a new file
rewrite_pdf(newData);
% rename
fn=ls('rw_*');
fn=fn(4:end-1);
eval(['!mv rw_',fn,' hb_',fn])

%% BIU users
% better clean accelerometers first if you have them 
% use createCleanFile with 'HeartBeat',0, like this
fileName = 'c,rfhp0.1Hz';
p=pdf4D(fileName);
cleanCoefs = createCleanFile(p, fileName,'byLF',256 ,'Method','Adaptive','xClean',[4,5,6],'byFFT',0,'HeartBeat',0);
% then you get a file called 'xc,lf_c,rfhp0.1Hz' or something similar
% then clean the heartbeat like this (it lookx for xc,lf_* file)
newData=correctHB; % or newData=correctHB([],[],1); if you want lots of plots
rewrite_pdf(newData);



