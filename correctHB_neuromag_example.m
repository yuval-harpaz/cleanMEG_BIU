% Example of cleaning neuromag data of HeartBeat artifact
% I read the fif file with fieldtrip
% If you don't want to use FieldTrip or if you want to rewrite a fif file
% see correctLF_neuromag_example.m

LS=ls('*.fif');
LS=LS(1:end-1);
hdr=ft_read_header(LS);
trl=[1,hdr.nSamples,0];
cfg=[];
cfg.trl=trl;
cfg.demean='yes';
cfg.dataset=LS;
cfg.channel='MEGMAG';
mag=ft_preprocessing(cfg);
meanMAG=mean(mag.trial{1,1});
cfg.channel='MEG';
meg=ft_preprocessing(cfg);
figOptions.label=meg.label;
figOptions.layout='neuromag306mag.lay';
clear mag
[cleanMEG,tempMEG,periodMEG,mcgMEG,RtopoMEG]=correctHB(meg.trial{1,1},meg.fsample,figOptions,meanMAG);
