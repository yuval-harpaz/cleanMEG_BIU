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
