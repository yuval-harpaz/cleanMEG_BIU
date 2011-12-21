function [F1,AllPSDmeg1,AllPSDmeg2]=ffter(pat1,BadChans,pat2)
% retorns the power density of a file or two and plots the mean MEG after
% ignoring bad channels

pdf1 = pdf4D(pat1);

% BadChans=[74;204];
%% pdf1 (dirty)

hdr = pdf1.header;
chim = channel_index(pdf1, 'meg', 'name');
chir = channel_index(pdf1, 'ref', 'name');

% lat = lat2ind(pdf1, 1, [0 hdr.epoch_data{1,1}.epoch_duration])
display('reading data')
MEG1=read_data_block(pdf1,[1 hdr.epoch_data{1,1}.pts_in_epoch],chim);
%REF1=read_data_block(pdf1,[1 hdr.epoch_data{1,1}.pts_in_epoch],chir);
samplingRate=get(pdf1,'dr');


display('ffting');
[AllPSDmeg1, F1] = allSpectra(MEG1,samplingRate,0.25,'FFT');

%[AllPSDref1, F1] = allSpectra(REF1,samplingRate,0.25,'FFT');
m=meanWnan(AllPSDmeg1,2); % finding bad channels
zeroChans=find(isnan(m));
if ~isempty(zeroChans);
    zeroChans %#ok<NOPRT>
end
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
if exist('pat2','var')
    pdf2 = pdf4D(pat2);
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
end
end