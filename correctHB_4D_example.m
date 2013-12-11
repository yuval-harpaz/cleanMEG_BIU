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
% use createCleanFile with 'HeartBeat',0,
% then you get a file called 'xc,lf_c,rfhp0.1Hz' or something similar
% then clean the heartbeat like this
fn='xc,lf_c,rfhp0.1Hz';
p=pdf4D(fn);
sRate=double(get(p,'dr'));
hdr = get(p, 'header');
nSamp=hdr.epoch_data{1,1}.pts_in_epoch;
chi = channel_index(p, 'meg', 'name');
data = read_data_block(p,[1 nSamp],chi);
data=double(data);
label=channel_label(p,chi)';
figOptions.label=label;
figOptions.layout='4D248.lay';
[cleanData,temp2e,period2,MCG,Rtopo]=correctHB(data,sRate,figOptions);
% or without the figOptions if you don't like them all over the screen


