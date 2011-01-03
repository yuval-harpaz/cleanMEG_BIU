function showSignal(mMEG, samplingRate, figH, QRS)
% show the cleaned mMEG

% Dec-2010 MA

figure(figH)

%% filter the mMEG
BandPassSpecObj=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
2,4,50,100,60,1,60,samplingRate);
BandPassFilt=design(BandPassSpecObj  ,'butter');
mMEGf = myFilt(mMEG,BandPassFilt);

%% convolve with QRS
range = round(0.05*samplingRate); % 50 ms
sMEG = conv(QRS,mMEGf);
sMEG = sMEG(range:end-range+1);

%% plot
plot(sMEG,'r');
legend('Original', 'HB cleaned')

return
