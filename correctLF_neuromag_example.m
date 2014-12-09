% Example of cleaning neuromag data of HeartBeat artifact
% I read the fif file with fiff_read_raw_segment, not in this repo
% I use 'time' chanLF because there is no ref chan and finding the zero
% crossing on mag channels or eog was not as good (Boston data).

infile='sub01_raw.fif';
outfile='sub01LF_raw.fif';
raw=fiff_setup_read_raw(infile);
[data,~] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp,1:306);
lastChans= fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp,307:312);
cfg=[];
cfg.Ncycles=1000;
clean=correctLF(data,raw.info.sfreq,'time',cfg);
clear data
saveas(1,'sub01LF.tiff');
copyfile(infile,outfile);
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
fiff_write_raw_buffer(outfid,[clean;lastChans],cals);
fiff_finish_writing_raw(outfid);
