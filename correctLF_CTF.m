%cfg.dataset='0_MEG000_20140909_Empty600Hz.meg4';
%% empty room
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
cfgLF.method='adaptive';
clean=correctLF(meg.trial{1,1},meg.fsample,ref.trial{1,1}(noisei,:),cfgLF);
close;
% test PSD only after some 20sec
samp20s=round(20*meg.fsample);
figure;
plot(mean(abs(fftBasic(meg.trial{1,1}(:,samp20s:end),meg.fsample))),'r');
hold on
plot(mean(abs(fftBasic(clean(:,samp20s:end),meg.fsample))),'g');
legend('mean raw','mean clean')
