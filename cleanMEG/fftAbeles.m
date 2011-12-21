pdf1 = pdf4D('c1,rfhp0.1Hz');
%pdf2 = pdf4D('c,rfDC,lp');
pdf2 = pdf4D('xc,lf_c1,rfhp0.1Hz');
%pdf3 = pdf4D('/opt/msw/data/megmas_data0/A0noise/c,noW6x/26.10.10@09:35/3/c,rfhp1.0Hz');
BadChans=[74;204];
%% pdf1 (dirty)

hdr = pdf1.header;
chim = channel_index(pdf1, 'meg', 'name');
chir = channel_index(pdf1, 'ref', 'name');

% lat = lat2ind(pdf1, 1, [0 hdr.epoch_data{1,1}.epoch_duration])
MEG1=read_data_block(pdf1,[1 hdr.epoch_data{1,1}.pts_in_epoch],chim);
%REF1=read_data_block(pdf1,[1 hdr.epoch_data{1,1}.pts_in_epoch],chir);
samplingRate=get(pdf1,'dr');



[AllPSDmeg1, F1] = allSpectra(MEG1,samplingRate,0.25,'FFT');

%[AllPSDref1, F1] = allSpectra(REF1,samplingRate,0.25,'FFT');
m=meanWnan(AllPSDmeg1,2); % finding bad channels
if isempty(BadChans);BadChans=find(m>(6*median(meanWnan(AllPSDmeg1,2))));end;BadChans
psd1=AllPSDmeg1;
psd1(BadChans,:)=NaN; % excluding Bad Channles
figure

semilogx(F1,sqrt(meanWnan(psd1)),'r');

xlabel ('Frequency Hz')

ylabel('SQRT(PSD), T/sqrt(Hz)');

title('Mean PSD for all channels')

%axis tight
hold on

%% pdf2 (clean)
hdr = pdf2.header;
chim = channel_index(pdf2, 'meg', 'name');
chir = channel_index(pdf2, 'ref', 'name');

% lat = lat2ind(pdf2, 1, [0 hdr.epoch_data{1,1}.epoch_duration])
MEG2=read_data_block(pdf2,[1 hdr.epoch_data{1,1}.pts_in_epoch],chim);
%REF2=read_data_block(pdf2,[1 hdr.epoch_data{1,1}.pts_in_epoch],chir);
samplingRate=get(pdf2,'dr');



AllPSDmeg2 = allSpectra(MEG2,samplingRate,0.25,'FFT');

%[AllPSDref2, F2] = allSpectra(REF2,samplingRate,0.25,'FFT');
psd2=AllPSDmeg2;
psd2(BadChans,:)=NaN;
%figure;

% semilogx(F2,sqrt(meanWnan(psd2)));
shiftPsd2=meanWnan(psd2);
shiftPsd2(1,2:end)=shiftPsd2(1:(end-1)); %shifts the second (blue) plot 1 point to the right to prevent overlap of curves
semilogx(F1,sqrt(shiftPsd2));
%%
% pdf3 = pdf4D('/opt/msw/data/megmas_data0/A0noise/c,noW6x/26.10.10@09:35/5/c,rfhp1.0Hz');
% hdr = pdf3.header;
% chim = channel_index(pdf3, 'meg', 'name');
% %chir = channel_index(pdf3, 'ref', 'name');
% MEG3=read_data_block(pdf3,[1 hdr.epoch_data{1,1}.pts_in_epoch],chim);
% samplingRate=get(pdf3,'dr');
% AllPSDmeg3 = allSpectra(MEG3,samplingRate,0.25,'FFT');
% psd3=AllPSDmeg3;
% psd3(BadChans,:)=NaN;
% figure;
% semilogx(F1,sqrt(meanWnan(psd1)),'r');
% hold on
% xlabel ('Frequency Hz')
% ylabel('SQRT(PSD), T/sqrt(Hz)');
% title('Mean PSD for all channels')
% shiftPsd3=meanWnan(psd3);
% shiftPsd3(1,2:end)=shiftPsd3(1:(end-1)); %shifts the second (blue) plot 1 point to the right to prevent overlap of curves
% semilogx(F1,sqrt(shiftPsd3));