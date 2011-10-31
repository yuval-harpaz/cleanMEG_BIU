function cleanCoefs = createCleanFileF(Tilda, inFile, varargin)
% Create a new pdf file from an old one
% cleanCoefs = createCleanFileE(pdf, inFile, , varargin);
%    e.g.
% >> inFile = 'C:\DATA\MEG\Noise\28-Oct-2008\onBed\c_rfhp1.0Hz';
% >> pdf = pdf4D(inFile);
% >> [cleanCoefs, ok] = createCleanFileE(pdf, inFile, 'byLF , 'byFFT');
%
% pdf     - any pdf4D oblect - not truly used here
% inFile  - full path and name of the old file
%   VARARGIN allows for parameter pairs ('name',value) to be specified
%   Legal pairs are
% 'outFile' - full path and name of the new file, or
%           - [] when "automatic" file name is generated (this is the
%           default)
% 'byLF'    - either [], or 0, or powerof 2 numeric values.
%      if the value is: 0 do not clean the Line frequency artefacts
%           2^n the value of the trigger bit for line f marker 
%          (e.g. 256  or  [256, 1024])
%      if []: 'byLF' will search for a bit which is time locked within 1
%           cycle either to 50 or 60 Hz. (default 'byLF',[])
%           This will be applied to MEG and REF channels and if available also to
%           EEG and external channels, this bit will then be removed from
%           the trig channel.
%      if more then one value (e.g. [256,1024]), then clean first by the
%          first bit, then is there is an apreciable difference between the
%          two trig-periods use also the second.
% 'Method'   - relevant only if cleaning by 'FFT'.  3 possibiliities:
%      'Global' - use the avergae LF locked cycle averaged over the entire
%          piece of data [this is the default].
%    'Adaptive' - Update the average cycle with a time constant of ~173 
%          cycles. Start with the average over the first 250 cycles
%    'phasePrecession' - compute the phase of each trigger mark and use it
%          (Not recommended unless the interference is very precise in
%          period and amplitude.)
% 'byFFT'    - 0 Do not use REF channels to clean the  MEG
%              1 If to clean also the MEG channels to obtain minimal
%                       power at all frequencies (default 0)
% 'RInclude' - 'All'       - means consider all REF channels
%            - 'Automatic' - means exclude channels with small dinamic
%                            range, (this is the default)
%            - array       - list of channel numbers to include
% 'RiLimit'  - value in [0,1], setting the dynamic range of REF to be
%              accepted.  (lower values mean be MORE selective). 
%              [default 0.85].
% 'HeteroCoefs' - to use previously prepared coeficients for cleaning
%              either []  - for computing coeficients with every piece
%              or coefStruct As obtained by coefsAllByFFT. (default [])
% 'DefOverflow' - define the value for overflow in MEG (default 1e-11).
% 'CleanPartOnly' - pair of times [T1,T2] in seconds. leave some parts 
%              as is and use only the part between T1 and T2.
% 'HeartBeat' - 0 donot attempt to clean heart beats (this is the default)
%               [] clean heart beat and find its period from the data
%               any value - the expected heart beat period in seconds.
%                 cleaning the mean HB after all the other cleaning.  This
%                 applies for both MEG, REF, and EEG if recorded.
% 'xClean'     - 0 donot use external inputs for cleaning
%                [1,2,...] list of external channels to use.  [default 0]
% 'xBands'     - list of frequency bounderies for cleaning
%                [] means use [1:139,140:20:800];
% 'chans2ignore' - list of channel numbers to be ignored.
% 'Memory'     - either 'small' or 'large'. or some number (e.g.100000).
%                How many samples of data to be processed in one iteration
% % NOT to be used.  Use xClean instead
% % chans2ignore' - ignore these channels. The value is
% %         the list of channels, (e.g. [74, 204]) which are bad and should
% %         not be treated.
%                
%
% cleanCoefs  - cell Array with coefsAllByFFT for each section 
%               (containing section definition, PCs for REF and 
%               weights for cleaning as well as the the Heart Beat).
%
% noQuestions - 1 - the program will NOT ask for user input (will overide
%                   exisiting clean file)
%               0 - the default, user input might be reqiured
% 
% NOTE this will work properly only if data in the file is continuous.
%      Due to memory limitations the file is read and written several
%      times. Each time the processing is in pieces.  First time the times
%      and mean waveshape of the Line Frequency is computed.  In the second
%      run the LF artifacts are cleaned,the position of heart beats is noted
%      and the mean HB is computed for each channel.
%      in the third pass the mean HB is subtracted, extarnal channels are 
%      used to clean the MEG and the file is cleaned by FFT. 

% Nov-2009  MA  Generated from CreateCleanFile on 15-Nov-2009
%  Updates
% Dec-2009 Algorithm for cleaning Heart Beat improved
% Jan-2010 MEG is filtered 1-100 before extracting the heartBeat times
% Oct-2010 XTR channels used to clean acceleration, LF cleaning improved
% Nov-2010 If heartBeats are too small or irregular they are not cleaned
% Dec-2010 When huge overflow is found (>1e-8) the data from there on is
%          converted to 0. problems with high heart rate and big noise on 
%          bad channels -fixed MA
% Dec-2010 Yossi 
%          replaced ~ with Tilda 
%          use fileparts() instead of manual strtok & flip
%          use copyfile() instead of manual strings & commands manipulation
%          new optional parameter 'noQuestions'

%% Initialize
% setup some default values
hugeVal = 1e-8;    % impossibly large MEG data
Method = 'Global'; % use global mean for cleaning LF
doXclean = false;  % default - no cleaning by XTR
xBand=[];
% External = false;  % no cleaning by XTR channels
outLierMargin = 20;  % Channels exceeding this value are ignored
overFlowTh= 1e-9;
minVarAcc = 1e-4;    % variance of acceleration must be above that
if ispc % adjust this according to the max memory in the system
    samplesPerPiece = 120000;
else
    samplesPerPiece = 500000;
end

warning('ON', 'MEGanalysis:missingInput:settingBands')
warning('ON','MATLAB:MEGanalysis:NotUsed')

%BTi Data Formats:
SHORT   =	1;
LONG    =	2;
FLOAT   =	3;
DOUBLE  =	4;

%% test for VARARGIN
if nargin > 2
    if rem((nargin-2),2)~= 0
        error('MATLAB:MEGanalysis:WrongNumberArgs',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'outFile', 'byLF', 'byFFT', 'RInclude', 'RiLimit','HeteroCoefs',...
        'DefOverflow', 'CleanPartOnly', 'HeartBeat', 'xClean', 'xBands', ...
        'Method', 'chans2ignore', 'Memory', 'noQuestions'};
    okargs = lower(okargs);
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error('MATLAB:MEGanalysis:BadParameter',...
                'Unknown parameter name:  %s.',pname);
        elseif length(k)>1
            error('MATLAB:MEGanalysis:BadParameter',...
                'Ambiguous parameter name:  %s.',pname);
        else
            switch(k)
                case 1  % outFile
                    outFile = pval;
                case 2  % byLF
                    byLF = pval;
                case 3 %  byFFT
                    doFFT = pval;
                case 4 %  RInclude
                    if ischar(pval)
                        RIval={'all','automatic'};
                        kk = strmatch(lower(pval), RIval);
                        if isempty(kk)
                            error('MATLAB:MEGanalysis:BadParameter',...
                                'Unknown parameter value:  %s.',pval);
                        elseif length(k)>1
                            error('MATLAB:MEGanalysis:BadParameter',...
                                'Ambiguous parameter value:  %s.',pval);
                        else
                            switch(k)
                                case 1  % None
                                    Rchans2Include = [];  % do not ignore any
                                    makeIgnoreList = false;
                                case 2  % Automatic
                                    makeIgnoreList = true;
                            end  % end of switch for RIval
                        end  % end of test for unique value
                    else  %A numeric list
                        Rchans2Include = pval;
                        makeIgnoreList = false;
                    end  % end of test for char values
                    case 5 %  RiLimit
                        RiLimit = pval;
                case 6 % 'HeteroCoefs'
                    coefStruct = pval;
                case 7 % 'DefOverflow'
                    overFlowTh = pval;
                case 8  % 'CleanPartOnly'
                    tStrt = pval(1);
                    tEnd  = pval(2);
                    aPieceOnly=true;
                case 9 % 'HeartBeat'
                    HBperiod = pval;
                    % doHB=true;
                case 10  % xClean
                    if isempty(pval) || any((pval==0))
                        doXclean = false;
                        xBand = [];
                    else
                        doXclean = true;
                        xChannels = pval;
                    end
                case 11  % xBand
                    if doXclean
                        xBand = pval;
                        if isempty(xBand)  % do the default
                            xBand = [1:139,140:20:800]';
                        end
                    end
                case 12  % method for LF cleaning
                    Method = upper(pval);
                    legalArgs = {'Global','Adaptive','phasePrecession'};
                    legalArgs = upper(legalArgs);
                    whichArg = strmatch(Method, legalArgs);
                    if ~isempty(whichArg)
                        switch whichArg
                            case 1 
                                Global=true;
                                Adaptive = false;
                                phasePrecession = false;
                            case 2
                                Global=false;
                                Adaptive = true;
                                phasePrecession = false;
                            case 3
                                Global=false;
                                Adaptive = false;
                                phasePrecession = true;
                        end  % end of switch
                    else
                        error ('MATLAB:MEGanalysis:IllegalParam'...
                            ,['Allwed METHODs are: ' legalArgs])
                    end
                case 13  %chans2ignore
                    chans2ignore = pval;
                case 14  % how much memory to use per chunk
                    if strcmpi(pval,'SMALL'), samplesPerPiece = 120000;
                    elseif strcmpi(pval,'LARGE'), samplesPerPiece = 500000;
                    elseif isnumeric(pval), samplesPerPiece = pval;
                    else 
                        error('MATLAB:pdf4D:ImproperParam',...
                        'allowed values for memory are: ''Small'', ''Large''')
                    end
                case 15  %noQuestions
                    noQuestions = pval;
            end  % end of switch
        end  % end of tests for unique arg name
    end  % end of testing for even number of argumants
end  % end of more then two input arguments


%% define the missing variables
if doXclean && ~exist('xBand','var') % define the xBands
     xBand = [1:139,140:20:800]';
end
if doXclean && isempty(xBand) % define the xBands
     xBand = [1:139,140:20:800]';
end
QRSstartT = -0.042;  % duration of QRS around peak in HB artifact
QRSendT   =  0.040;
clipAmplitude = true; % do not let relative amplitudes go beyond limits
amplitudeFudge=0.8;
if ~exist('chans2ignore', 'var'), chans2ignore=[]; end
if ~exist('outFile','var'), outFile =[]; end
if ~exist('RiLimit','var'), RiLimit =[]; end
if isempty(RiLimit), RiLimit = 0.85; end 
if ~exist('overFlowTh', 'var'), overFlowTh=[]; end
if isempty(overFlowTh), overFlowTh= 1e-9; end
if ~exist('doFFT','var'), doFFT =[]; end
if isempty(doFFT) 
    doFFT = false;
end
if ~exist('byLF','var'), byLF =[]; end
if isempty(byLF)
    doLineF = true; 
    findLF = true;
elseif byLF==0
        findLF = false;
        doLineF = false;
else % not a char - assume a number
    doLineF = true;
    lineF = byLF;
    findLF = false;
end
if ~exist('makeIgnoreList','var'), makeIgnoreList =true; end
if ~exist('coefStruct','var'), coefStruct=[]; end
if isempty(coefStruct)
    externalCoeficients = false;
else
    externalCoeficients = true;
end
if ~exist('aPieceOnly','var'), aPieceOnly=false; end
if ~exist('tStrt','var'), tStrt=[]; end
if isempty(tStrt), tStrt=0; end
if ~exist('HBperiod','var'), HBperiod=0; end
if isempty(HBperiod)
    doHB=true;
    findHBperiod=true;
elseif HBperiod>0
    doHB=true;
    findHBperiod=false;
else
    doHB=false;
end
legalArgs = {'Global','Adaptive','phasePrecession'};
legalArgs = upper(legalArgs);
whichArg = strmatch(upper(Method), legalArgs);
if ~isempty(whichArg)
    switch whichArg
        case 1
            Global=true;
            Adaptive = false;
            phasePrecession = false;
        case 2
            Global=false;
            Adaptive = true;
            phasePrecession = false;
        case 3
            Global=false;
            Adaptive = false;
            phasePrecession = true;
    end  % end of switch
end
if ~exist('noQuestions', 'var')
    noQuestions = 0;
end
    
pIn=pdf4D(inFile);
%% get basic info on the file
hdr = get(pIn,'header');
samplingRate   = double(get(pIn,'dr'));
numEpochs = length(hdr.epoch_data);
% lastSample = double(hdr.epoch_data{1}.pts_in_epoch);
if numEpochs>1
    error ('MATLAB:MEGanalysis:ImproperData',...
        'For now createCleanFileE doesnot work with Epoched data, try createCleanFile.')
    % warning ('MATLAB:pdf4D:notContinuous','epoched data, data "stitched" at bounderies')
%     epoched = true;
%     epochStart= nan(1,numEpochs);
%     lastSample=0;
%     for ii=1:numEpochs
%         epochStart(ii) = lastSample+1;
%         lastSample = lastSample+double(hdr.epoch_data{ii}.pts_in_epoch);
%     end
%     epochEnds = [epochStart(2:end)-1, lastSample];
else
    lastSample = double(hdr.epoch_data{1}.pts_in_epoch);
    epochs = [1,lastSample];
end
% totalT = lastSample/samplingRate;

%% prepare for some cleaning
if ~exist('tEnd','var'), tEnd=[]; end
if isempty(tEnd), tEnd=lastSample/samplingRate; end
maxF = 0.4*samplingRate; % do not process beyond that
% dF = 1/samplingRate;
if doXclean  % truncate the top frequencies
    for ii=1:length(xBand)
        if xBand(ii) > maxF
            break
        end
    end
    xBand(ii)=maxF;
    xBand(ii+1:end)=[];
end
if numEpochs>1
    if doHB
        error('MATLAB:pdf4D:notContinuous',...
            'Cannot clean Heart Beat on epoched files')
    end
    if tStrt>0
        error ('MATLAB:pdf4D:notContinuous',...
            'Cleaning only part of data doesnot work on epoched recordings')
    end
end

%% read a piece of the trig signal to find if Line frequency is available
testT=30;  % 30 seconds
% linePeriod = round(samplingRate/50);  % the default value
if lastSample/samplingRate<testT
    testT = lastSample/samplingRate;
end
if doLineF 
    chit = channel_index(pIn,'TRIGGER');
    trig = read_data_block(pIn,double(samplingRate*[1,testT]),chit);
    if isempty(trig)
        warning('MATLAB:MEGanalysis:noData','Line Frequency trig not found')
        doLineF = false;
    elseif findLF
        lineF = findLFbit(trig,samplingRate);
        if isempty(lineF)
            warning('MATLAB:MEGanalysis:noData','Line Frequency trig not found')
            doLineF = false;
        else
            doLineF = true;
        end
    end
end
if doLineF
    whereUp=find(diff(mod(trig,2*lineF(1)-1)>=lineF(1))==1);
    if isempty(whereUp)
        doLineF=false;
        warning ('MEGanalysis:missingParam','Couldnot clean the line artefacts')
        linePeriod=10;
    else
        linePeriod = round(mean(diff(whereUp)));
    end
    if length(lineF)>1  % two trig bits are supplied
        whereUp2=find(diff(mod(trig,2*lineF(2)-1)>=lineF(2))==1);
        if isempty(whereUp2)
            warning ('MEGanalysis:missingParam',...
                '2nd bit for LF is not flipping -- ignored')
        else
            linePeriod2 = mean(diff(whereUp2));
            if abs(linePeriod2-mean(diff(whereUp)))> 0.5/samplingRate
                useTrig2=true;
                linePeriod2=round(linePeriod2);
            else
                useTrig2=false; 
                lineF = lineF(1);
            end
        end
    else
        useTrig2 = false;
        linePeriod2=[];
    end
    df=1/linePeriod;
    whereUpAll = nan(1,ceil(1.2*(lastSample/linePeriod)));
    whereUpNext=1;

end

%% prepare list of MEG channels to ignore
chi = channel_index(pIn,'meg');
chn = channel_name(pIn,chi);
numMEGchans = length(chi);
ignoreChans = false(1, numMEGchans);
for ii = 1:length(chans2ignore)
    ignoreChans(chans2ignore(ii)) = true;
end
% make a list of valid channels
[Tilda, chiSorted] = sortMEGnames(chn,chi);
validChans = chiSorted;
validChans(chans2ignore) = [];

%% HB period

if doHB
    mMEG = mean(read_data_block(pIn,double(samplingRate*[1,testT]),validChans),1);
%     % define the filter
%     BP4to50Spec=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
%           2,4,50,100,60,1,60,samplingRate);
%     BP4to50Filt=design(BP5to50Spec  ,'butter');
%     mMEGf=myFilt(mMEG, BP4to50Filt);
    [whereisHB, zTime, Errors]= findHB01(mMEG, samplingRate,HBperiod,...
        'PLOT', 'VERIFY');
    if (sum(Errors.shortHB)>3) || (sum(Errors.longHB)>3) || (Errors.numSmallPeaks>3) % added by Yuval
        warning('MATLAB:MEGanalysis:ImproperData',...
            'Heart beat signal not clear enough No HB cleaning done');
        disp(Errors)
        doHB = false;
    end % end of added by Yuval
    if findHBperiod  % find the period
        HBperiod = 0.5*(max(diff(whereisHB))+min(diff(whereisHB)))/samplingRate;
    end
    HBcycleLength = round(HBperiod*samplingRate);
    % from here on zTime, iBefore and iAfter are FIXED!
    iBefore = zTime;
    iAfter = HBcycleLength -iBefore-1;    
end

%% find which REF channels to ignore
chirf = channel_index(pIn,'ref');
if makeIgnoreList
    %  algorithm
    % for each channel:
    %    make a histogram of REF values with 256 bins
    %    find the 10% to 90% values
    %    if the fraction of bins with 0 counts there is < RiLimit accept
    %    this channel.
    REF = read_data_block(pIn,double(samplingRate*[1,testT]),chirf);
    % search for channels with low dinamic range
    fN=zeros(1,length(chirf));
    ilst = 1:length(chirf);
    for jj = 1:length(ilst)
        ii = ilst(jj);
        hst=hist(diff(REF(ii,:)),256);
        chst = cumsum(hst);
        lo=find(chst>0.1*chst(end),1);
        hi=find(chst>0.9*chst(end),1);
        if lo==hi
            f=1;
        else
            f = sum(hst(lo:hi)==0)/(hi-lo);
        end
        fN(jj)=f;
    end
    Rchans2Include = find(fN<=RiLimit);
    if length(Rchans2Include)<9
        sFN = sort(fN);
        newRL = sFN(9);
        Rchans2Include = find(fN<=newRL);
        warning('MATLAB:MEGanalysis:notEnoughInputs', ...
            ['Not enough REF channels with enough dynamic range!' ...'
            '\n changing threshold to %d'],newRL)
        % doFFT=false;
    end
    numREFchans = length(Rchans2Include);
end
%% prepare the output structure
cleanInfo = struct('fileName',inFile , 'samplingRate',samplingRate ,...
    'allCoefs',[] , 'eigVec',[] , 'bands',[], 'xClean',doXclean ,...
    'xBand',xBand);
% prepare the list of indices for heart beat
if doHB
    HBparams = struct('whereisHB',[] , 'Errors',[] , 'zTime',[] , 'mMEGhb',[] ...
        ,'MEGamplitude',[], 'iBefore',[] , 'iAfter',[] , 'firstT',[] ...
        ,'lastT',[] , 'REFamplitude',[] ...
        ,'rChannels',Rchans2Include , 'mREFhb',[]);
    
end

%% continue checkups
if ~doLineF && ~doFFT && ~doHB && ~doXclean % && ~External
    warning('MATLAB:MEGanalysis:notEnoughInputs', 'Nothing to clean - ABORTING!')
    return
end
if isempty(outFile)
    whichF=0;
    % break into Path and Name
    %[Name,Path] = strtok(fliplr(inFile),'\'); %% this is a bug fix by
    %Yossi
%     Name = fliplr(Name);
%     Path = fliplr(Path);
    [Path, name, ext] = fileparts(inFile);
    Name = [name ext];
    
    if doHB,  whichF = whichF+1; end
    if doFFT, whichF = whichF+2; end
    if doLineF, whichF = whichF+4; end
    if doXclean, whichF = whichF+8; end
    % if External, whichF = whichF+16; end
    switch whichF
        case 1
            outFpreFix ='hb';
        case 2
            outFpreFix ='cf';
        case 3
            outFpreFix = 'cf,hb';
        case 4
            outFpreFix = 'lf';
        case 5
            outFpreFix = 'hb,lf';
        case 6
           outFpreFix = 'cf,lf';
        case 7
            outFpreFix = 'cf,hb,lf';
        case 8
            outFpreFix = 'xc';
        case 9
            outFpreFix = 'xc,lf';
        case 10
            outFpreFix = 'xc,cf';
        case 11
            outFpreFix = 'xc,cf,hb';
        case 12
            outFpreFix = 'xc,lf';
        case 13
            outFpreFix = 'xc,hb,lf';
        case 14
            outFpreFix = 'xc,cf,lf';
        case 15
            outFpreFix = 'xc,cf,hb,lf';
%         case 16
%             outFpreFix = 'xt';
%         case 17
%             outFpreFix = 'xt,hb';
%         case 18
%             outFpreFix ='xt,cf';
%         case 19
%             outFpreFix = 'xt,cf,hb';
%         case 20
%             outFpreFix = 'xt,lf';
%         case 21
%             outFpreFix = 'xt,hb,lf';
%         case 22
%            outFpreFix = 'xt,cf,lf';
%         case 23
%             outFpreFix = 'xt,cf,hb,lf';
%         case 24
%             outFpreFix = 'xt,xc';
%         case 25
%             outFpreFix = 'xt,xc,lf';
%         case 26
%             outFpreFix = 'xt,xc,cf';
%         case 27
%             outFpreFix = 'xt,xc,cf,hb';
%         case 28
%             outFpreFix = 'xt,xc,lf';
%         case 29
%             outFpreFix = 'xt,xc,hb,lf';
%         case 30
%             outFpreFix = 'xt,xc,cf,lf';
%         case 31
%             outFpreFix = 'xt,xc,cf,hb,lf';
            
    end
    % generate a new file name
    [S,R] = strtok(Name,'_');  % find if there is already an underscore
    if isempty(R)
        Name = [outFpreFix, '_' S];
    else
       Name = [outFpreFix, ',' Name];
    end
    outFile = [Path Name];
end
if ~noQuestions
    if exist(outFile, 'file')
        R = input(['outFile ' outFile ' already exists. \n Do you wish to replace? [1-yes/0-no] ']);
        if R~=1
            disp('Aborting')
            return
        end
    end
end

copyfile(inFile, outFile);
% if ispc
%     command = ['copy "' inFile '" "' outFile '"'];
%     disp(['Wait while : ' command])
%     [s,w]=dos(command);
%     if s~=0  % error
%         error('MATLAB:MEGanalysis:fileNames',['could not copy files: ' w])
%     end
% elseif isunix
%     command = ['cp ' inFile ' ' outFile];
%     unix(command);
% else
%     warning('MATLAB:MEGanalysis:unknownSystem','Unsupported operating system - nothing done')
%     return
% end

%% start the file and get parameters
disp('Starting the MEG file and adjusting parameters')
p=pdf4D(outFile);
% cnf = p.config;
hdr = get(p,'header');
%empty header means no pdf
if isempty(hdr)
    error('MATLAB:MEGanalysis:noPDFfile','Need pdf to write data')
end

numEpochs = length(hdr.epoch_data);
if numEpochs>1
    error ('MATLAB:MEGanalysis:ImproperData',...
        'For now createCleanFileE doesnot work with Epoched data, try createCleanFile.')
%     % error('MATLAB:pdf4D:notContinuous','Cannot clean epoched files')
%     warning ('MATLAB:pdf4D:notContinuous','epoched data, data "stitched" at bounderies')
%     epoched = true;
%     epochStart= nan(1,numEpochs);
%     lastSample=0;
%     for ii=1:numEpochs
%         epochStart(ii) = lastSample+1;
%         lastSample = lastSample+double(hdr.epoch_data{ii}.pts_in_epoch);
%     end
%     epochEnds = [epochStart(2:end)-1, lastSample];
else
    epoched = false;
    lastSample = double(hdr.epoch_data{1}.pts_in_epoch);
end
samplingRate   = double(get(p,'dr'));
chi = channel_index(p,'meg');
chn = channel_name(p,chi);
% numMEGchans = length(chi);
[Tilda, chiSorted] = sortMEGnames(chn,chi);
chix = channel_index(p,'EXTERNAL');
chit = channel_index(p,'TRIGGER');
% chir = channel_index(p,'RESPONSE');
if ~isempty(chix)
    chnx = channel_name(p,chix);
    [~, chixSorted] = sortMEGnames(chnx,chix);
end
chie = channel_index(p,'EEG');
if ~isempty(chie)
    che = channel_name(p,chie);
    [Tilda, chieSorted] = sortMEGnames(che,chie);
    numEEGchans = length(che);
end
chirf = channel_index(p,'ref');
if ~isempty(Rchans2Include)
    chirf = chirf(Rchans2Include);
    % numREFchans = length(chirf);
end

%% readjust the bounderies
if samplesPerPiece>lastSample
    samplesPerPiece= lastSample;
end
if ~epoched
    if  ~aPieceOnly
        numPieces = ceil(lastSample/samplesPerPiece);
        samplesPerPiece= floor(lastSample/numPieces);
        startApiece = 1:samplesPerPiece:lastSample;
        stopApiece  = startApiece+samplesPerPiece;
        deltaEnd = lastSample-stopApiece(end);
        if deltaEnd<0
            stopApiece(end) = lastSample;
        else
            error ('wrong division of data')
        end
    else
        numSamples = round((tEnd-tStrt)*samplingRate);
        numPieces = ceil(numSamples/samplesPerPiece);
        samplesPerPiece= floor(numSamples/numPieces);
        firstS = floor(tStrt*samplingRate);
        lastSample = firstS + numSamples-1;
        startApiece = firstS:samplesPerPiece:lastSample;
        stopApiece  = startApiece+samplesPerPiece;
        deltaEnd = lastSample-stopApiece(end);
        if deltaEnd<0
            stopApiece(end) = lastSample;
        else
            error ('wrong division of data')
        end
    end
end

%% decide on size of time slice to process
% find if a line frequency trigger exists
if ~ doLineF
    df = 1/round(0.02*samplingRate);
end
transitionFactors = (0.5*df:df:1-0.5*df); % factors for merging near the cutoff between files
numTransition = length(transitionFactors);
transitionFactors = repmat(transitionFactors,numMEGchans,1);

%% prepare for reading and writing the file
%total number of channels in pdf
total_chans = double(hdr.header_data.total_chans);

%open file (to read and write), always big endean
fid = fopen(outFile, 'r+', 'b');

if fid == -1
    error('Cannot open file %s', outFile);
end

%% clean the data
%  From here on we pass over the entire file up to 4 times.
% 1.  If doLF we go over the entire file without writing back anything and
%     mark the times of LF trig + compute the mean LF cycle over the entire
%     file.
% 2.  If doLF subtract the mean LF cycle.  
%     If doHB find the time of all heart-beats and the mean HB shape. If doLF
%     write the LF-cleaned data back to the file
% 4.  If two LF trigs and they have different periods, repeat (1) again for
%     the 2-nd trig over the entire file, then read the file subtract it
%     and write back to file.
% 3.  If doHB go over the entire file, clean the mean HB,
%     If doXchannel clean by FFT the external channels from the MEG
%     If doFFT clean the data by REFs. and save back.

if doLineF
    LFcycleLength = max(diff(whereUp)+3);
    % meanPeriod = mean(whereUp);
    if Global || Adaptive  % use one average for the entire file
        MEGlfCycle = zeros(numMEGchans,LFcycleLength);
        REFlfCycle = zeros(length(chirf),LFcycleLength);
        if exist('chieSorted','var')
            EEGlfCycle = zeros(length(chieSorted),LFcycleLength);
        end
        if exist('chixSorted','var')
            XTRlfCycle = zeros(length(chixSorted),LFcycleLength);
        end
    elseif PhasePrecession % use interpolation
        MEGlfCycle = zeros(numMEGchans,3*LFcycleLength);
        REFlfCycle = zeros(length(chirf),3*LFcycleLength);
        if exist('chieSorted','var')
            EEGlfCycle = zeros(length(chieSorted),3*LFcycleLength);
        end
        if exist('chixSorted','var')
            XTRlfCycle = zeros(length(chixSorted),3*LFcycleLength);
        end
    else
        disp(['The only legal cleaning methods are: ''Global'''...
            ', ''Adaptive'', or ''PhasePrecession'''])
    end
    totalCycles =0; % cycles of LF
end

%% initialize the structure array of HB
if doHB
%     for ii = 1:numPieces
        % HBparams.firstT = startApiece(ii)/samplingRate;
    HB(1:numPieces) = HBparams;
%     end
    minAmplitude=1;
    maxAmplitude=1;
end
if doFFT
    cleanCoefs(1:numPieces) = cleanInfo;
else
    cleanCoefs=[];
end
    
%%  Pass 1 find LF cycles and mean cycle
if doLineF
    for ii = 1:numPieces
        startI = startApiece(ii);
        endI  = stopApiece(ii);
        disp(['finding LF for: ' num2str([startI,endI]/samplingRate)])
        
        trig = read_data_block(p, [startI,endI], chit);
        whereUp=find(diff(mod(trig,2*lineF-1)>=lineF)==1);
        
        % read all types of data
        MEG = read_data_block(p, [startI,endI], chiSorted);
        % ignoreChans = false(1,size(MEG,1));
        chans2analyze = find(~ignoreChans);
        for iii=chans2analyze % 1:size(MEG,1)
            if sum(MEG(iii,:)>1e-9)>100
                warning('MEGanalysis:overflow',...
                    ['In channel #' num2str(iii) ' more then 100 oveflows - ignored'])
                % MEG(iii,:)=0;
                ignoreChans(iii)=true;
            end
        end
        % check if any channel has too big variance
        vv = var(MEG,[],2);
        outVv = (vv-median(vv))/mad(vv);
        mVv = max(outVv);
        if mVv > outLierMargin
            ovIndx = find(outVv>outLierMargin);
            for jj = 1:length(ovIndx)
                iii = ovIndx(jj);
                if ~ignoreChans(iii)  % add to ignored list
                    warning('MEGanalysis:tooBig',...
                        ['In channel #' num2str(iii) ' Variance is too big - ignored'])
                    % MEG(iii,:)=0;
                    ignoreChans(iii)=true;
                end
            end
        end
        % sometimes the last value is HUGE??
        [Tilda,junkData] = find(MEG(~ignoreChans,:)>hugeVal,1);
        if ~isempty(junkData)
            endI = startI + size(MEG,2)-1;
            warning('MATLAB:MEGanalysis:nonValidData', ...
                ['MEGanalysis:overflow','Some MEG values are huge at: ',...
                num2str(endI/samplingRate) ' - truncated'])
            MEG(:,junkData:end)=0;
            % timeList(end) = junkData-1;
            % else
            % timeList(end) = size(MEG,2);
        end
        
        if exist('chieSorted','var')
            EEG = read_data_block(p, [startI,endI], chieSorted);
        end
        if ~exist('chirf','var')
            chirf = channel_index(p,'ref');
        end
        if ~isempty(chirf)  % read the reference channels
            % chnrf = channel_name(p,chirf);
            REF = read_data_block(p, [startI,endI], chirf);
        end
        if exist('chixSorted','var')  % read the reference channels
            XTR = read_data_block(p, [startI,endI], chixSorted);
        end
        % Check if clear by acceleration was required
        if doXclean
            if any(var(XTR(xChannels,:),[],2)<minVarAcc)
                error('MATLAB:MEGanalysis:ImproperParam',...
                    'Not enough data in external channels for cleaning')
            end
        end
        % find where are the LF cycles
        numCyclesHere = length(whereUp);
        whereUpAll(whereUpNext:whereUpNext+numCyclesHere-1) = whereUp+startI-1;
        whereUpNext = whereUpNext +numCyclesHere;
        if Global % find the overall mean LF cycle
            maxL = LFcycleLength-1;
            sumC = zeros(numMEGchans,maxL+1);
            for cycle = 1:numCyclesHere-1;
                startCycle = whereUp(cycle);
                if startCycle+maxL <= size(MEG,2)
                    sumC = sumC + MEG(:,startCycle:startCycle+maxL);
                end
            end
            MEGlfCycle = MEGlfCycle+sumC;
            totalCycles = totalCycles+numCyclesHere-1;
            sumC = zeros(numREFchans,maxL+1);
            for cycle = 1:numCyclesHere-1;
                startCycle = whereUp(cycle);
                if startCycle+maxL <= size(MEG,2)
                    sumC = sumC + REF(:,startCycle:startCycle+maxL);
                end
            end
            REFlfCycle = REFlfCycle+sumC;
            if exist('chieSorted','var')
                sumC = zeros(size(EEG,1),maxL+1);
                for cycle = 1:numCyclesHere-1;
                    startCycle = whereUp(cycle);
                    if startCycle+maxL <= size(MEG,2)
                        sumC = sumC + EEG(:,startCycle:startCycle+maxL);
                    end
                end
                EEGlfCycle = EEGlfCycle+sumC;
            end
            if exist('chixSorted','var')
                sumC = zeros(size(XTR,1),maxL+1);
                for cycle = 1:numCyclesHere-1;
                    startCycle = whereUp(cycle);
                    if startCycle+maxL <= size(MEG,2)
                        sumC = sumC + XTR(:,startCycle:startCycle+maxL);
                    end
                end
                XTRlfCycle = XTRlfCycle+sumC;
            end
        end % end of finding global LF cycles
        clear MEG REF startI endI trig
        if exist('chieSorted','var'),  clear EEG; end
        if exist('chixSorted','var'),  clear XTR; end
        if doLineF, clear trig whereUp; end
    end  % end of going over all pieces
    
    if Global
        MEGlfCycle = MEGlfCycle/totalCycles;
        MEGlfCycle = MEGlfCycle - repmat(mean(MEGlfCycle,2),1,maxL+1);
        REFlfCycle = REFlfCycle/totalCycles;
        REFlfCycle = REFlfCycle - repmat(mean(REFlfCycle,2),1,maxL+1);
        if exist('chieSorted','var')
            EEGlfCycle = EEGlfCycle/totalCycles;
            EEGlfCycle = EEGlfCycle - repmat(mean(EEGlfCycle,2),1,maxL+1);
        end
        if exist('chixSorted','var')
            XTRlfCycle = XTRlfCycle/totalCycles;
            XTRlfCycle = XTRlfCycle - repmat(mean(XTRlfCycle,2),1,maxL+1);
        end
    end % end of normalizing global averages
end

%% pass 2 go again over the data, clean and save back
transitions=nan(1,numPieces);
if doHB||doLineF||doXclean
    if Adaptive % oreoare arrays for ongoing mean LF cycles
        if doLineF
            meanL = max(diff(whereUpAll)) +1;
            meanMEG = zeros (numMEGchans,meanL);
            meanREF = zeros (numREFchans,meanL);
            if exist('chieSorted','var')
                meanEEG = zeros (numEEGchans,meanL);
            end
             if exist('chixSorted','var')
                 meanXTR = zeros (size(chixSorted,2),meanL);
             end
        end
    end
    for ii = 1:numPieces
        startI = startApiece(ii);
        endI  = stopApiece(ii);
        if doLineF
            disp(['cleaning LF for the piece ' num2str([startI,endI]/samplingRate)])
        end
        transitions(ii) = endI;
        MEG = read_data_block(p, [startI,endI], chiSorted);
        % sometimes the last value is HUGE??
        [Tilda,junkData] = find(MEG>hugeVal,1);
        if ~isempty(junkData)
            endI = startI + size(MEG,2)-1;
            warning('MATLAB:MEGanalysis:nonValidData', ...
                ['MEGanalysis:overflow','Some MEG values are huge at: ',...
                num2str(endI/samplingRate) ' - truncated'])
            MEG(:,junkData:end)=0;
        end
        
        if exist('chieSorted','var')
            EEG = read_data_block(p, [startI,endI], chieSorted);
        end
        if ~exist('chirf','var')
            chirf = channel_index(p,'ref');
        end
        if ~isempty(chirf)  % read the reference channels
            % chnrf = channel_name(p,chirf);
            REF = read_data_block(p, [startI,endI], chirf);
        end
        if exist('chixSorted','var')  % read the reference channels
            XTR = read_data_block(p, [startI,endI], chixSorted);
        end
        
        %%    % clean the 50Hz if needed
        if doLineF
            whereUp = whereUpAll(whereUpAll>startI & whereUpAll<endI)-startI+1;
            for nc = chans2analyze
                x=MEG(nc,:);
                unCleaned=[];
                if Global
                    [y,unCleaned] = cleanMean(x, whereUp, MEGlfCycle(nc,:), 1, []);
                elseif Adaptive
                    if ii ==1 % first time
                        [y, mean1] = cleanLineF(x, whereUp,...
                            epochs, 'Adaptive');
                        % define the array of means
                        meanMEG(nc,:) = mean1;
                    else
                        mean1 = meanMEG(nc,:);
                        [y, mean1] = cleanLineF(x, whereUp,...
                            epochs, 'Adaptive', mean1);
                        meanMEG(nc,:) = mean1;
                    end
                elseif phasePrecession % must be phase precession
                    y = cleanLineF(x, whereUp, epochs, 'phasePrecession');
                else
                   error ('MATLAB:MEGanalysis:IllegalParam'...
                            ,['Allwed METHODs are: ' legalArgs])
                end
                MEG(nc,:)=y;
            end
            if ~isempty(unCleaned)
                unCleaned = (strtI+unCleaned)/samplingRate;
                for uc = 1:length(unCleaned)
                    disp(['Could not clean a HB cycle at ', num2str(unCleaned(uc))]);
                end
            end
            for nc = 1:length(chirf)  % not all REF chans are used
                x=REF(nc,:);
                if Global
                    y = cleanMean(x, whereUp, REFlfCycle(nc,:), 1, []);
                elseif Adaptive
                    if ii ==1 % first time
                        [y, mean1] = cleanLineF(x, whereUp,...
                            epochs, 'Adaptive');
                        % define the array of means
                        meanREF(nc,:) = mean1;
                    else
                        mean1 = meanREF(nc,:);
                        [y, mean1] = cleanLineF(x, whereUp,...
                            epochs, 'Adaptive', mean1);
                        meanREF(nc,:) = mean1;
                    end
                elseif phasePrecession % must be phase precession
                    y = cleanLineF(x, whereUp, epochs, 'phasePrecession');
                else
                   error ('MATLAB:MEGanalysis:IllegalParam'...
                            ,['Allwed METHODs are: ' legalArgs])
                end
                REF(nc,:)=y;
            end
            if exist('chieSorted','var')  % clean the 50 Hz
                for nc = 1:numEEGchans
                    x=EEG(nc,:);
                    if Global
                        y = cleanMean(x, whereUp, EEGlfCycle(nc,:), 1, []);
                    elseif Adaptive
                        if ii ==1 % first time
                            [y, mean1] = cleanLineF(x, whereUp,...
                                epochs, 'Adaptive');
                            % define the array of means
                            meanEEG(nc,:) = mean1;
                        else
                            mean1 = meanEEG(nc,:);
                            [y, mean1] = cleanLineF(x, whereUp,...
                                epochs, 'Adaptive', mean1);
                            meanEEG(nc,:) = mean1;
                        end
                    elseif phasePrecession % must be phase precession
                        y = cleanLineF(x, whereUp, epochs, 'phasePrecession');
                    else
                        error ('MATLAB:MEGanalysis:IllegalParam'...
                            ,['Allwed METHODs are: ' legalArgs])
                    end
                    EEG(nc,:)=y;
                end
            end
            if exist('chixSorted','var')  % clean the 50 Hz
                for nc = 1:size(chixSorted,2)
                    x=XTR(nc,:);
                    if Global
                        y = cleanMean(x, whereUp, XTRlfCycle(nc,:), 1, []);
                    elseif Adaptive
                        if ii ==1 % first time
                            [y, mean1] = cleanLineF(x, whereUp,...
                                epochs, 'Adaptive');
                        % define the array of means
                        meanXTR(nc,:) = mean1;
                        else
                            mean1 = meanXTR(nc,:);
                            [y, mean1] = cleanLineF(x, whereUp,...
                                epochs, 'Adaptive', mean1);
                            meanXTR(nc,:) = mean1;
                        end
                    elseif phasePrecession % must be phase precession
                        y = cleanLineF(x, whereUp, epochs, 'phasePrecession');
                    else
                        error ('MATLAB:MEGanalysis:IllegalParam'...
                            ,['Allwed METHODs are: ' legalArgs])
                    end
                    XTR(nc,:)=y;
                end
            end
        end
        if useTrig2 || ~isempty(linePeriod2)
            warning('MATLAB:MEGanalysis:NotUsed','Second LF trig is not used for now.')
            warning('OFF','MATLAB:MEGanalysis:NotUsed');  % do not repeat this warning
        end
        %% clean by the extarnal lines if needed
        if doXclean
            disp(['cleaning XTR channels for the piece ' num2str([startI,endI]/samplingRate)])
            for nc = chans2analyze
                x=MEG(nc,:);
                y = coefsAllByFFT(x, XTR(xChannels,:),...
                    samplingRate, xBand, [], 1,false);
                MEG(nc,:) = y;
            end
            for nc = 1:numREFchans
                x=REF(nc,:);
                y = coefsAllByFFT(x, XTR(xChannels,:),...
                    samplingRate, xBand, [], 1,false);
                REF(nc,:) = y;
            end
        end
        %% search for heart beats
        if doHB
            disp(['Finding heart beat for the piece ' num2str([startI,endI]/samplingRate)])
            mMEG = mean(MEG(~ignoreChans,:),1);
%             mMEGf=myFilt(mMEG, BP5to50Filt);
            [whereisHB, Tilda, Errors, amplitudes] = findHB01(mMEG, samplingRate, HBperiod,'NOPLOT','Verify');
            amplitudes = amplitudes/mean(amplitudes);
            if sum(amplitudes<=0)  > 0
                error('MATLAB:MEGanalysis:improperValue',...
                    'Negative HB amplitudes at: %d to %d', startI, endI)
            end
            if min(amplitudes)<minAmplitude
                minAmplitude = min(amplitudes);
            end
            if max(amplitudes)>maxAmplitude
                maxAmplitude = max(amplitudes);
            end
            % find the mean for channel
            mMEGhb = zeros(numMEGchans,HBcycleLength);
            for chan = 1:numMEGchans
                x=MEG(chan,:);
                mMEGhb(chan,:) = meanAround(x, whereisHB,iBefore, iAfter);
            end
            HBparams.whereisHB = whereisHB;
            HBparams.Errors = Errors;
            HBparams.zTime = zTime;
            HBparams.mMEGhb = mMEGhb;
            HBparams.iBefore = iBefore;
            HBparams.iAfter = iAfter;
            HBparams.firstT = startI/samplingRate;
            % repeat for the REFs
            mREFhb = zeros(numREFchans,HBcycleLength);
            for chan = 1:numREFchans
                x=REF(chan,:);
                mREFhb(chan,:) = meanAround(x, whereisHB,iBefore, iAfter);
            end
            HBparams.mREFhb = mREFhb;
            % store
            HB(ii) = HBparams;
        end  % end of finding Heart beat positions
        
        if doLineF||doXclean
            %% write the cleaned from 50Hz and XTR data back
            %% replace the old data by the cleaned data
            switch hdr.header_data.data_format
                case SHORT
                    data_format = 'int16=>int16';
                    data_format_out = 'int16';
                    time_slice_size = 2 * total_chans;
                    config = get(p, 'config');
                    if isempty(config)
                        error('No Config: Could not scale data\n');
                    end
                    scale = channel_scale(config, hdr, 1:total_chans);
                    MEG = int16(MEG ./ repmat(scale(chiSorted)', 1, size(MEG,2)));
                    REF = int16(REF ./ repmat(scale(chirf)', 1, size(REF,2)));
                    if exist('chieSorted','var')  % EEG exists
                        EEG = int16(EEG ./ repmat(scale(chie)', 1, size(EEG,2)));
                    end
                    if exist('chixSorted','var')  % XTR exists
                        XTR = int16(XTR ./ repmat(scale(chix)', 1, size(XTR,2)));
                    end
                case LONG
                    data_format = 'int32=>int32';
                    data_format_out = 'int32';
                    time_slice_size = 4 * total_chans;
                    config = get(p, 'config');
                    if isempty(config)
                        error('No Config: Could not scale data\n');
                    end
                    scale = channel_scale(config, hdr, 1:total_chans);
                    MEG = int32(MEG ./ repmat(scale(chi)', 1, size(MEG,2)));
                    REF = int32(REF ./ repmat(scale(chirf)', 1, size(REF,2)));
                    if exist('chieSorted','var')  % EEG exists
                        EEG = int32(EEG ./ repmat(scale(chie)', 1, size(EEG,2)));
                    end
                    if exist('chixSorted','var')  % XTR exists
                        XTR = int32(XTR ./ repmat(scale(chix)', 1, size(XTR,2)));
                    end
                case FLOAT
                    data_format = 'float32=>float32';
                    data_format_out = 'float32';
                    time_slice_size = 4 * total_chans;
                    MEG = single(MEG);
                    REF = single(REF);
                    if exist('chieSorted','var')  % EEG exists
                        EEG = single(EEG);
                    end
                    if exist('chixSorted','var')  % XTR exists
                        XTR = single(XTR);
                    end
                case DOUBLE
                    data_format = 'double';
                    time_slice_size = 8 * total_chans;
                    MEG = double(MEG);
                    REF = double(REF);
                    if exist('chieSorted','var')  % EEG exists
                        EEG = double(EEG);
                    end
                    if exist('chixSorted','var')  % XTR exists
                        XTR = double(XTR);
                    end
                otherwise
                    error('Wrong data format : %d\n', hdr.header_data.data_format);
            end
            %skip some time slices
            lat = startI;
            %     if lat>1  %  skip
            status = fseek(fid, time_slice_size * (lat-1), 'bof');
            if status~=0
                error('MEGanalysis:pdf:fileOperation', ['Did not advance the file ' ferror(fid)])
            end
            %     end
            % Read the old data and replace with new
            oldData = fread(fid, [total_chans, size(MEG,2)], data_format);
            oldData(chiSorted,:) = MEG;  % replace the MEG channels
            oldData(chirf,:) = REF;  % replace the REF channels
            if exist('chieSorted','var')
                oldData(chieSorted,:) = EEG;  % replace the EEG channels
            end
            if exist('chixSorted','var')
                oldData(chixSorted,:) = XTR;  % replace the XTR channels
            end
            if doLineF
                trig=oldData(chit,:);
                trig = clearBits(trig, lineF);
                oldData(chit,:)=trig;
            end
            status = fseek(fid, time_slice_size * (lat-1), 'bof');
            if status~=0
                error('MEGanalysis:pdf:fileOperation', ['Did not advance the file ' ferror(fid)])
            end
            
            fwrite(fid, oldData, data_format_out);
        end
    end

    %% clean the space for next group
    clear MEG REF oldData startI endI trig

    if exist('chieSorted','var'),  clear EEG; end
    if exist('chixSorted','var'),  clear XTR; end
    if doLineF, clear trig whereUp; end
    if doHB, clear x; end
end  % end of treating one piece for LF 

if ~doFFT && ~doHB
    return
end

if doHB
    minAmplitude = amplitudeFudge*minAmplitude;
    maxAmplitude = maxAmplitude/amplitudeFudge;
    % extract the means
    % to avoid inclusion of part of QRS at the end of a cycle we limit the mean
    % beat length to the mid between and longest sortest Hb cycle over all 
    % the data (provided it is not too short)
    % we also find what is the maximal no. of HB per section of data for
    % preparing the Amplitudes arrays
    minHBcycle = size(mREFhb,2);
    maxHBcycle = minHBcycle;
    maxNoHB=0;
    for ii = 1:length(HB)
        minHere = min(diff(HB(ii).whereisHB));
        maxHere = max(diff(HB(ii).whereisHB));
        if maxHere>maxHBcycle
            maxHBcycle = maxHere;
        end
        if minHere > 1.8*zTime % use it
            % minHBcycle = minHere;
        else
            whereShort = whereisHB(diff(HB(ii).whereisHB)==minHere)/samplingRate;
            warning('MATLAB:MEGanalysis:inapropriateValue','Shortest cycle too small at %d',...
                whereShort);
            % add to the errors
            err = HB(ii).Errors;
            shErr = [err.shortHB whereShort];
            HB(ii).Errors.shortHB = shErr;
        end
        if length(HB(ii).whereisHB)>maxNoHB
            maxNoHB = length(HB(ii).whereisHB);
        end
    end
    % HBcycle = HBcycleLength;  %round((minHBcycle+maxHBcycle)/2);
%     iBefore = zTime;
%     iAfter = HBcycle-zTime-1;
    cycleLength = iBefore+iAfter+1;
    MEGhbCycle = zeros(size(mMEGhb,1),cycleLength);
    REFhbCycle = zeros(size(mREFhb,1),cycleLength);
    numCycles = 0;
    for ii=1:length(HB);
        numHere = length(HB(ii).whereisHB);
        rangeHere = 1:cycleLength;
        numCycles = numCycles +numHere;
        MEGhbCycle = MEGhbCycle +numHere*HB(ii).mMEGhb(:,rangeHere);
        REFhbCycle = REFhbCycle +numHere*HB(ii).mREFhb(:,rangeHere);
    end
    MEGhbCycle = MEGhbCycle/numCycles;
    REFhbCycle = REFhbCycle/numCycles;
    % offset so edges are near zero
    for jj = 1:size(MEGhbCycle,1)
        tailOffset = mean([MEGhbCycle(jj,1:10) MEGhbCycle(jj,end-9:end)]);
        MEGhbCycle(jj,:) = MEGhbCycle(jj,:) - tailOffset;
    end
    for jj = 1:size(REFhbCycle,1)
        tailOffset = mean([REFhbCycle(jj,1:10) REFhbCycle(jj,end-9:end)]);
        REFhbCycle(jj,:) = REFhbCycle(jj,:) - tailOffset;
    end
end

%% clean HB and by FFT
for ii = 1:numPieces
    startI = startApiece(ii);
    endI  = stopApiece(ii);
    if doHB&&doFFT
        disp(['cleaning HB and by FFT for the piece ' num2str([startI,endI]/samplingRate)])
    elseif doHB &&~doFFT
        disp(['cleaning HB for the piece ' num2str([startI,endI]/samplingRate)])
    else
        disp(['cleaning by FFT for the piece ' num2str([startI,endI]/samplingRate)])
    end
    transitions(ii) = endI;

    % read all types of data
    MEG = read_data_block(p, [startI,endI], chiSorted);
    % sometimes the last value is HUGE??
    [Tilda,junkData] = find(MEG>hugeVal,1);
    if ~isempty(junkData)
        endI = startI + size(MEG,2)-1;
        warning('MATLAB:MEGanalysis:nonValidData', ...
            ['MEGanalysis:overflow','Some MEG values are huge at: ',...
            num2str(endI/samplingRate) ' - truncated'])
        MEG(:,junkData:end)=0;
    end
    
    if ~exist('chirf','var')
        chirf = channel_index(p,'ref');
    end
    if ~isempty(chirf)  % read the reference channels
        % chnrf = channel_name(p,chirf);
        REF = read_data_block(p, [startI,endI], chirf);
    end
    if exist('chieSorted','var')  % EEG exists
        EEG = read_data_block(p, [startI,endI], chieSorted);
    end
    if exist('chixSorted','var')  % read the reference channels
        XTR = read_data_block(p, [startI,endI], chixSorted);
    end
    if doHB
        % whereisHB has to be recomputed because the start of each section
        % is shifted so the start of susequent pieces overlap
        allAmplMEG = zeros(numMEGchans, maxNoHB);
        allAmplREF = zeros(numREFchans, maxNoHB);
        mMEG = mean(MEG(chans2analyze,:),1);
%         mMEGf=myFilt(mMEG, BP4to50Filt);
        whereisHB = findHB01(mMEG, samplingRate,[], 'noPLOT', 'noVERIFY');
        % check for a missing beat at the edges
        if whereisHB(1)>max(diff(whereisHB))  % missed at start
            tmpHB = [whereisHB(1)-mean(diff(whereisHB)); whereisHB];
            whereisHB = tmpHB;
        end
        if length(mMEG)-whereisHB(end)>max(diff(whereisHB))  % missed at end
            tmpHB = [whereisHB; whereisHB(end)+mean(diff(whereisHB))];
            whereisHB = tmpHB;
        end
        QRSbefore = -round(QRSstartT*samplingRate);
        QRSafter = round(QRSendT*samplingRate);
        % extract the parameters of HB for this piece of data
        HB(ii).whereisHB = whereisHB;
        % Errors = HB(ii).Errors;
        HB(ii).zTime = zTime;
        % mMEGhb = HB(ii).mMEGhb;
        HB(ii).firstT = startI/samplingRate;
        HB(ii).lastT = endI/samplingRate;
        % mREFhb = HB(ii).mREFhb;
        QRSstartI = zTime - QRSbefore;
        QRSendI   = zTime + QRSafter;
        % for each MEG channel find the Amplitudes of fit to the mean QRS
        % of that channel
        numChans = size(mMEGhb,1);
        numBeats = length(whereisHB);
        % cycleLength = size(mMEGhb,2);
        % Use only the large QRS shape to normalize amplitude.
        nTmplt = MEGhbCycle(:, QRSstartI:QRSendI);
        nTmplt = nTmplt-repmat(mean(nTmplt,2),1,size(nTmplt,2));
        for jj=1:numChans
            x= nTmplt(jj,:);
            nTmplt(jj,:) = x/sqrt(x*x');
        end
        nTmplt = nTmplt'; % each normalized template in one column
        
        % test if first/lst not too close to th edge
        if whereisHB(1)-QRSbefore <=0
            jStrt=2;
        else
            jStrt=1;
        end
        if whereisHB(end)+QRSafter >size(MEG,2)
            jEnd=numBeats-1;
        else
            jEnd=numBeats;
        end
        
        for chan = 1:numChans
            Amplitudes = zeros(size(whereisHB));
            x=MEG(chan,:);
            for jj = jStrt:jEnd  %BUT take care of first and last cycles!!
                % what if too near the begining? =XXX
                thisBeat = x((whereisHB(jj)-QRSbefore):(whereisHB(jj)+QRSafter));
                thisBeat = thisBeat-mean(thisBeat);
                Amplitudes(jj) = thisBeat*nTmplt(:,chan);
            end
            Amplitudes(jStrt:jEnd) = Amplitudes(jStrt:jEnd)/mean(Amplitudes(jStrt:jEnd));
            if jStrt>1
                Amplitudes(1)=1;
            end
            if jEnd>numBeats
                Amplitudes(end)=1;
            end
            % clip amplitudes
            if clipAmplitude
                Amplitudes(Amplitudes<minAmplitude) = minAmplitude;
                Amplitudes(Amplitudes>maxAmplitude) = maxAmplitude;
            end
            allAmplMEG(chan,1:length(Amplitudes)) = Amplitudes;
            y = cleanMean(x, whereisHB, MEGhbCycle(chan,:), zTime, Amplitudes);
            MEG(chan,:)=y;
        end
        HB(ii).MEGamplitude = allAmplMEG;
        % clean also the REF channels
        nTmplt = REFhbCycle(:, QRSstartI:QRSendI);
        nTmplt = nTmplt-repmat(mean(nTmplt,2),1,size(nTmplt,2));
        for jj=1:numREFchans
            x= nTmplt(jj,:);
            nTmplt(jj,:) = x/sqrt(x*x');
        end
        nTmplt = nTmplt'; % each normalized template in one column
        for nc = 1:numREFchans
            Amplitudes = zeros(size(whereisHB));
            x=REF(nc,:);
            for jj = jStrt:jEnd  %BUT take care of first and last cycles!!
                thisBeat = x((whereisHB(jj)-QRSbefore):(whereisHB(jj)+QRSafter));
                thisBeat = thisBeat-mean(thisBeat);
                Amplitudes(jj) = thisBeat*nTmplt(:,nc);
            end
            Amplitudes(jStrt:jEnd) = Amplitudes(jStrt:jEnd)/mean(Amplitudes(jStrt:jEnd));
            if jStrt>1
                Amplitudes(1)=1;
            end
            if jEnd>numBeats
                Amplitudes(end)=1;
            end
            % clip amplitudes
            if clipAmplitude
                Amplitudes(Amplitudes<minAmplitude) = minAmplitude;
                Amplitudes(Amplitudes>maxAmplitude) = maxAmplitude;
            end
            allAmplREF(nc,1:length(Amplitudes)) = Amplitudes;
            y = cleanMean(x, whereisHB, REFhbCycle(nc,:), zTime, Amplitudes);
            REF(nc,:)=y;
        end
        HB(ii).REFamplitude = allAmplREF;
    end

    if doFFT
        if ~externalCoeficients
            sumOverflow = sum(sum(abs(MEG(chans2analyze,:))>overFlowTh))>0;
            if ii==1  % test for overflow on first piece
                if ~sumOverflow
                    [MEG(chans2analyze,:), coefStruct]  = ...
                        coefsAllByFFT(MEG(chans2analyze,:),REF,samplingRate,...
                        [],[],overFlowTh);
                else % cannot overcome the overflow problem
                    if ispc 
                        delete(outFile);
                    elseif isunix
                        command = ['rm' outFile];
                        unix(command);
                    else
                        disp(['Unknown system. Remove ' outFile ' manually!'])
                    end
                    error('MATLAB:MEGanalysis:cannotPerformFunction','cannot Clean. Cleaned file is removed.')
                end
            else  % ii>1 first section already cleaned, and auto-clean
                if ~sumOverflow
                    [MEG(chans2analyze,:), coefStruct]  = coefsAllByFFT(MEG(chans2analyze,:),REF,samplingRate,...
                        [],[], overFlowTh);
                else  % overflow & ii>1
                    [MEG(chans2analyze,:), coefStruct]  = coefsAllByFFT(MEG(chans2analyze,:),REF,samplingRate,...
                        coefStruct.bands, coefStruct, overFlowTh);
                end
            end  % end for test if ii==1
        else  % use external coeficients
            [MEG, coefStruct]  = coefsAllByFFT(MEG,REF,samplingRate,...
                coefStruct.bands, coefStruct, overFlowTh);
        end
        if ii==1 % the start (and end) of ifft may be large values
            upAmplitude = (0:10)/10;
            MEG(:,1:length(upAmplitude)) = MEG(:,1:length(upAmplitude))...
                     .*repmat(upAmplitude,size(MEG,1),1);
        end
        if ii>1  % make a smooth transition between the end of last piece
            % and the begining of the new one
            MEG(:,1:numTransition) = endOfPiece.*(1-transitionFactors) +...
                transitionFactors.*MEG(:,1:numTransition);
        end
        % save the last samples for creating a smooth transition to the
        % next piece
        if ii<numPieces % truncate the last piece
            endOfPiece = MEG(:,end-numTransition+1:end);
            MEG(:,end-numTransition+1:end) = [];
            REF(:,end-numTransition+1:end) = [];
%             trig(end-numTransition+1:end)=[];
            transitions(ii) = startI + size(MEG,2)-1;
            startApiece(ii+1) = transitions(ii)+1;  % Start next piece from 
            %                        the beginning of the overlapping piece
            if exist('chieSorted','var')  % EEG exists
                EEG(:,end-numTransition+1:end) = [];
            end
            if exist('chixSorted','var')  % XTR exists
                XTR(:,end-numTransition+1:end) = [];
            end
        end
        if ii==numPieces % the start (and end) of ifft may be large values
            dwnAmplitude = (10:-1:0)/10;
            numP = length(dwnAmplitude);
            MEG(chans2analyze,end-numP+1:end) = MEG(chans2analyze,end-numP+1:end)...
                .*repmat(upAmplitude,length(chans2analyze),1);
        end
    end

%% replace the old data by the cleaned data
    switch hdr.header_data.data_format
        case SHORT
            data_format = 'int16=>int16';
            data_format_out = 'int16';
            time_slice_size = 2 * total_chans;
            config = get(p, 'config');
            if isempty(config)
                error('No Config: Could not scale data\n');
            end
            scale = channel_scale(config, hdr, 1:total_chans);
            MEG = int16(MEG ./ repmat(scale(chiSorted)', 1, size(MEG,2)));
            REF = int16(REF ./ repmat(scale(chirf)', 1, size(REF,2)));
            if exist('chieSorted','var')  % EEG exists
                EEG = int16(EEG ./ repmat(scale(chie)', 1, size(EEG,2)));
            end
            if exist('chixSorted','var')  % XTR exists
                XTR = int16(XTR ./ repmat(scale(chix)', 1, size(XTR,2)));
            end
        case LONG
            data_format = 'int32=>int32';
            data_format_out = 'int32';
            time_slice_size = 4 * total_chans;
            config = get(p, 'config');
            if isempty(config)
                error('No Config: Could not scale data\n');
            end
            scale = channel_scale(config, hdr, 1:total_chans);
            MEG = int32(MEG ./ repmat(scale(chi)', 1, size(MEG,2)));
            REF = int32(REF ./ repmat(scale(chirf)', 1, size(REF,2)));
            if exist('chieSorted','var')  % EEG exists
                EEG = int32(EEG ./ repmat(scale(chie)', 1, size(EEG,2)));
            end
            if exist('chixSorted','var')  % XTR exists
                XTR = int32(XTR ./ repmat(scale(chix)', 1, size(XTR,2)));
            end
        case FLOAT
            data_format = 'float32=>float32';
            data_format_out = 'float32';
            time_slice_size = 4 * total_chans;
            MEG = single(MEG);
            REF = single(REF);
            if exist('chieSorted','var')  % EEG exists
                EEG = single(EEG);
            end
            if exist('chixSorted','var')  % XTR exists
                XTR = single(XTR);
            end
        case DOUBLE
            data_format = 'double';
            time_slice_size = 8 * total_chans;
            MEG = double(MEG);
            REF = double(REF);
            if exist('chieSorted','var')  % EEG exists
                EEG = double(EEG);
            end
            if exist('chixSorted','var')  % XTR exists
                XTR = double(XTR);
            end
        otherwise
            error('Wrong data format : %d\n', hdr.header_data.data_format);
    end
    %skip some time slices
    lat = startI;
    %     if lat>1  %  skip
    status = fseek(fid, time_slice_size * (lat-1), 'bof');
    if status~=0
        error('MEGanalysis:pdf:fileOperation', ['Did not advance the file ' ferror(fid)])
    end
    %     end
    % Read the old data and replace with new
    oldData = fread(fid, [total_chans, size(MEG,2)], data_format);
    oldData(chiSorted,:) = MEG;  % replace the MEG channels
    oldData(chirf,:) = REF;  % replace the REF channels
    if exist('chieSorted','var')
        oldData(chieSorted,:) = EEG;  % replace the EEG channels
    end
    if exist('chixSorted','var')
        oldData(chixSorted,:) = XTR;  % replace the XTR channels
    end
    status = fseek(fid, time_slice_size * (lat-1), 'bof');
    if status~=0
        error('MEGanalysis:pdf:fileOperation', ['Did not advance the file ' ferror(fid)])
    end

    fwrite(fid, oldData, data_format_out);
    
    %% clean the space for next group
    clear MEG REF oldData startI endI

    if exist('chieSorted','var'),  clear EEG; end
    if exist('chixSorted','var'),  clear XTR; end
%     if doLineF, clear whereUp; end
    if doHB
        HBparams = HB(ii);
        cleanCoefs(ii).HBparams = HBparams;
        clear x  whereisHB  Amplitudes mMEG AllAmplMEG AllAmplREF
        clear HBparams
    else
        cleanCoefs(ii).HBparams = [];
    end
    if doFFT
        cleanCoefs(ii).coefStruct = coefStruct;
        % clear coefStruct  % they should stay on in case of overflow
    else
        cleanCoefs(ii).coefStruct = [];
    end
end  % end of treating one piece for HB and FFT

%% wrap up

fclose(fid);
warning('ON', 'MEGanalysis:missingInput:settingBands')

return
