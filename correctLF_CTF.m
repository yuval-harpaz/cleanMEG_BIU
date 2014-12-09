% reading the MEG 
cfg=[];
cfg.dataset='0_MEG000_20140909_4kHz.ds/0_MEG000_20140909_4kHz.meg4';
cfg.channel = {'MEG'};
cfg.continuous = 'yes';
cfg.demean='yes';
meg=ft_preprocessing(cfg);
% reading REF channels
cfg.channel = {'MEGREF'};
ref=ft_preprocessing(cfg);
f=abs(fftBasic(ref.trial{1,1},ref.fsample));
[~,noisei]=max(f(:,60)); % choosing the noisiest ref channel for zero crossing detection
cfgLF=[];
cfgLF.Lfreq=60;
cfgLF.Ncycle=1000;
clean=correctLF(meg.trial{1,1},meg.fsample,ref.trial{1,1}(noisei,:),cfgLF);
