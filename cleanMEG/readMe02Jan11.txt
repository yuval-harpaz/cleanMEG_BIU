updates for createCleanFileH
19-07-2011
1.  The following argument was added:
    'hugeVal'  It controls when a channel is considered as unusable.
2.  A bug when runing the procedure from the same directory as the 
    data was fixed
3.  When there were several errors of the same type in the heart beat detection 
    it caused an ERROR message and the procedure stopped - fixed.
4.  When multiple bits are specified for LF cleaning the procedure tests whether 
    they fall within the 50 or 60 Hz range.  Those that are not are ignored.  
    Then if timing of two bits are synchronized within +/-1 sample only the 
    first one in the list is used.  (I hope in the future to introduce there 
    also a video-sync bit.)  Only if after all these two bits are left - a 
    Warning message is issued and the first bit alone is used.
    If the argument is empty ('byLF',[]), all the bits that are marking 50 or 60
    Hz are first detected.

updates for createCleanFileG
24-02-2011
1.  A few bugs in HB cleaning that caused it to fail occasionaly (when HB was 
    close to the end or the begining) - fixed.
2.  Big step like artifacts are detected and corrected for computing the mean HB.
3.  Big step like artifacts may be overcome in the data by a new option:
    'stepCorrect',1  The position of the step is notified. The correction is not 
    perfect a small "spike" remains at the start of the step.
4.  Only one trig bit is used for 50 Hz cleaning.  If you specify 2 bits the 2nd 
    one is ignored and a warning message is issued.
5.  Excessive trig - bits may be erased in the cleaned file. new option: 
    'maskTrigBits'
6.  If HB was cleaned the cleaned mean trace is plotted on the top of the initial 
    plot.
7.  If the file is big, the cleaning is done one piece at a time.  The size of this
    piece was increased.  If it causes "out of memory" error you may specify smaller 
    pieces by the option:
    'Memory'

I have tried the program on all the data pieces I had from Yossi and Odelia and 
it worked on all of them.

If you are runnig an older version of MATLAB and the symbol ~ is not known, replace 
all places where this is used as a name of variable by the word  Tilda. (But note 
that some of ~ is used to negate logic values - do not change these.)

Please try it carefully before using regularly
Moshe
