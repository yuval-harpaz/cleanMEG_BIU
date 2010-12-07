fileName = 'c1,rfhp0.1Hz';
p=pdf4D(fileName);
cleanCoefs = createCleanFileF(p, fileName,'byLF',256 ,'Method','Adaptive' , 'xClean',[4,5,6], 'chans2ignore',[74,204] ,'byFFT',0 , 'HeartBeat',0);