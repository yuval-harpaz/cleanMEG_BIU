function correctXC(fileName)
% designed to use createCleanFile for external channels (4D data) after
% correctLF and before correctHB
if ~exist('fileName','var')
    try
        LS=ls('hb*');
        fileName=LS(1:end-1)
    end
    if ~exist('LS','var')
        try
            LS=ls('lf*');
            fileName=LS(1:end-1)
        end
    end
    if ~exist('LS','var')
        try
            LS=ls('c,rf*');
            fileName=LS(1:end-1)
        catch
            error('cant get file name')
        end
    end
end
p=pdf4D(fileName);
cleanCoefs = createCleanFile(p, fileName,'byLF',0 ,'xClean',[4,5,6],'byFFT',0,'HeartBeat',0,'outLierMargin',40);
f0=median(abs(fftRaw(fileName)));
figure;
plot(f0(1:120),'r');
ylabel('PSD')
xlabel('Hz')
hold on
[Path, name, ext] = fileparts(fileName);
try
    f1=median(abs(fftRaw(fullfile(Path,['xc,',name,ext]))));
catch
    f1=median(abs(fftRaw(fullfile(Path,['xc_',name,ext]))));
end
%f1=median(abs(fftRaw(['xc,',fileName])));
plot(f1(1:120),'g');

