function compareHBclean(varargin)

% verify the effectivness of the heart beat cleaning by averaging to the
% onset of the heart-beats 
% use as:
% compareHBclean('-argument_name', argument_value)
% Inputs:
% '-file1'  - name of the uncleaned MEG file. Default = 'c,rfhp0.1Hz'.
% '-file2'  - name of the cleaned MEG file as given by createCleanFile.m.
%           Default = 'hb_c,rfhp0.1Hz'.
% '-file3'   - [] (default - no ICA comparison) or name of the MEG file for heart beat
%           cleaning using ICA.
% '-acc'    - check the average of the accelerometrs time locked to the
%           heart beats to check for balistocardiogram artifact.
%           [] - do not perform (default)
%           'reg' - regular average.
%           'abs' - average of the absolute value.
%           'both' - both regular average and average of the absolute
%           value.
% '-numTrials'  - number of trials to compare. Default = 250.
% '-accOnly'    - 1 - perform only accelerometers average without heart beat
%               cleaning comparison.
%               0 - Perform also heart beat comparison (Default).
% '-runICA'     - 1 - run ICA on the 3rd dataset
%               - 0 - load dataica from previous (Default)
% '-saveAvgs    - 1 - save avg1, avg2 and avg3
%               - 0 - do not save (Default)
% '-saveFigs'   - 1 - save the averages figures.
%               - 0 do not save (Default)
% '-ICA'        - 1 - perform comparisson with ICA
%               - 0 - do not perform comparisson with ICA (Default)

% NOTES:        - The fieldtrip toolbox and the function process_varargin.m should
%                be on your matlab path.
%               - the current folder should contain cleanCoefs.mat 
%               - When using with ica dataica.mat should be in the current
%               directory

%   IT Nov 2012


% process input arguments
okargs=        {'-file1'       '-file2'         '-file3'             '-acc'     '-numTrials' '-onlyAcc' '-runICA' '-saveAvgs' '-saveFigs', '-ICA'};
default_values={'c,rfhp0.1Hz', 'hb_c,rfhp0.1Hz', [],                 'none',     250,          0,        0,        0,          0,           0};
okargs = lower(okargs);
keyvalue = process_varargin('compareHBclean', okargs, default_values, varargin);
file1 = keyvalue{1};
file2 = keyvalue{2};
file3 = keyvalue{3};
acc   = keyvalue{4};
numTrials = keyvalue{5}+1;
accOnly   = keyvalue{6};
runICA = keyvalue{7};
saveAvgs = keyvalue{8};
saveFigs = keyvalue{9};
ICA = keyvalue{10};

% read the MEG file
p = pdf4D(file1);
% get the sampling rate
samplingRate   = double(get(p,'dr')); %get the sqampling rate
% read the cleanCoefs.mat as given by createCleanFile.m after heart beat
% cleaning
if exist('cleanCoefs.mat')
    load cleanCoefs
else
    error('cleanCoefs.mat file does not exist')
    return
end
% if the nuber of trials exceeds the number of trials for the first section
% take the number of trials for the first section only
if length(cleanCoefs(1).HBparams.whereisHB) < numTrials
    numTrials = length(cleanCoefs(1).HBparams.whereisHB);
end
% take the exact times of the heart beats
trlTrig = cleanCoefs(1).HBparams.whereisHB(2:numTrials);
%adjust number of trials
numTrials = length(trlTrig);


%% check ERP for accelerometers
if acc=='none';
    acc = [];
end
if ~isempty(acc)
    % define offset for trl matrix
    offset = round(-0.28*samplingRate);
    % define trial length for trl matrix
    trlLength = round(1*samplingRate);
    % define trl matrix
    trl = zeros(numTrials,3);
    trl(:,1) = trlTrig + offset;
    trl(:,2) = trl(:,1) + trlLength;
    trl(:,3) = offset;
    % read the accelerometers data using fieldtrip
    cfg = [];
    cfg.dataset = file1;
    cfg.trl = trl;
    cfg.channel ={'X4' 'X5' 'X6'};
    dataAcc = ft_preprocessing(cfg);
    
    if strcmp(acc, 'reg') % regular average
        % perform timelock analysis
        cfg = [];
        avgAcc = ft_timelockanalysis(cfg,dataAcc);
        % plot the results (butterfly plot)
        tMEG = avgAcc.time;
        figure
        plot(tMEG, avgAcc.avg(1,:),'b')
        hold on
        plot(tMEG, avgAcc.avg(2,:),'r')
        plot(tMEG, avgAcc.avg(3,:),'g')
        set(gca,'fontsize',14);
        title('accelerometers locked to HB')
        xlabel('Time (sec.)')
        ylabel('Amplitude');
        legend('X', 'Y', 'Z')
        grid on
        axis tight
    elseif strcmp(acc, 'abs') % absolute value
        
        % take the absolut value for all trials
        for ii = 1 : length(dataAcc.trial)
            dataAcc.trial{1,ii} = abs(dataAcc.trial{1,ii});
        end
        % perform timelock analysis
        cfg = [];
        avgAcc = ft_timelockanalysis(cfg,dataAcc);
        % plot the results
        tMEG = avgAcc.time;
        figure
        plot(tMEG, avgAcc.avg(1,:),'b')
        hold on
        plot(tMEG, avgAcc.avg(2,:),'r')
        plot(tMEG, avgAcc.avg(3,:),'g')
        set(gca,'fontsize',14);
        title('abs accelerometers locked to HB')
        xlabel('Time (sec.)')
        ylabel('Amplitude');
        legend('X', 'Y', 'Z')
        grid on
        axis tight
    elseif strcmp(acc, 'both') % take both average and absolute value
        
        % take the absolut value for all trials
        dataAccAbs = dataAcc;
        for ii = 1 : length(dataAcc.trial)
            dataAccAbs.trial{1,ii} = abs(dataAcc.trial{1,ii});
        end
        % calculate regular average
        cfg = [];
        avgAcc = ft_timelockanalysis(cfg,dataAcc);
        % calculate average of the absolute value
        cfg = [];
        avgAccAbs = ft_timelockanalysis(cfg,dataAccAbs);
        
        % plot the results of the regular average
        tMEG = avgAcc.time;
        figure
        plot(tMEG, avgAcc.avg(1,:),'b')
        hold on
        plot(tMEG, avgAcc.avg(2,:),'r')
        plot(tMEG, avgAcc.avg(3,:),'g')
        set(gca,'fontsize',14);
        title('accelerometers locked to HB')
        xlabel('Time (sec.)')
        ylabel('Amplitude');
        legend('X', 'Y', 'Z')
        grid on
        axis tight
        
        % plot the results of the absolute value
        tMEG = avgAcc.time;
        figure
        plot(tMEG, avgAcc.avg(1,:),'b')
        hold on
        plot(tMEG, avgAccAbs.avg(2,:),'r')
        plot(tMEG, avgAccAbs.avg(3,:),'g')
        set(gca,'fontsize',14);
        title('abs accelerometers locked to HB')
        xlabel('Time (sec.)')
        ylabel('Amplitude');
        legend('X', 'Y', 'Z')
        grid on
        axis tight
    end
    if accOnly % if heart beat comparison is not to be performed
        return
    end
end

%% start HB comparisson
% define offset for trl matrix
offset = round(-0.3*samplingRate);
% define trial length for trl matrix
trlLength = round(1*samplingRate);
% define trl matrix
trl = zeros(numTrials,3);
trl(:,1) = trlTrig + offset;
trl(:,2) = trl(:,1) + trlLength;
trl(:,3) = offset;
%% first dataset - before cleaning
% read data using fieldtrip
cfg = [];
cfg.dataset = file1;
cfg.trl = trl;
cfg.channel ={'MEG'};
cfg.demean = 'yes';
data1 = ft_preprocessing(cfg);
% calculate average
cfg = [];
cfg.channel = {'MEG'};
avg1 = ft_timelockanalysis(cfg,data1);

%% second dataset - after cleaning
% read the data using fieldtrip
cfg = [];
cfg.dataset = file2;
cfg.trl = trl;
cfg.channel = {'MEG'};
cfg.demean = 'yes';
data2 = ft_preprocessing(cfg);
% calculate average
cfg = [];
cfg.channel = {'MEG'};
avg2 = ft_timelockanalysis(cfg,data2);


%% Third dataset - Run ICA
if isempty(file3) % if not performing ICA comparison
    % plot the rsults of heart beat cleaning with createCleanFile.m
    tMEG = avg1.time;
    figure
    plot(tMEG, avg1.avg(1,:),'b')
    set(gca,'fontsize',14);
    hold on
    plot(tMEG, avg2.avg(1,:),'r')
    legend('None', 'Current Method');
    plot(tMEG, avg1.avg,'b')
    hold on
    plot(tMEG, avg2.avg,'r')
    
    title('compare HB clean all channels')
    xlabel('Time (sec.)')
    ylabel('Amplitude (Tesla)');
    grid on
    
elseif ICA % perform ICA comparison
    
    fprintf('========== Loading dataica from previous run ============');
    
    load dataica
    
    % calculate average timelocked to th heart beats positions
    cfg = [];
    cfg.channel = {'MEG'};
    avg3 = ft_timelockanalysis(cfg,dataica);
    
    % plot the results for ICA comparison
    tMEG = avg1.time;
    figure
    plot(tMEG, avg1.avg(1,:),'b')
    set(gca,'fontsize',14);
    hold on
    plot(tMEG, avg2.avg(1,:),'g')
    plot(tMEG, avg3.avg(1,:),'r')
    legend('avg1', 'avg2','avg3');
    plot(tMEG, avg1.avg,'b')
    hold on
    plot(tMEG, avg2.avg,'g')
    plot(tMEG, avg3.avg,'r')
    
    title('compare HB clean all channels')
    xlabel('Time (sec.)')
    ylabel('Amplitude (Tesla)');
    grid on
    if saveAvgs
        
        save avgs avg1 avg2 avg3
        
    end
    if saveFigs
        
        saveas(gcf, 'butterAll.fig');
        saveas(gcf, 'butterAll.png');
        
    end
else  % no ICA but the file is not empty - compare different cleaning
    
    % read the data using fieldtrip
    cfg = [];
    cfg.dataset = file3;
    cfg.trl = trl;
    cfg.channel = {'MEG'};
    cfg.demean = 'yes';
    data3 = ft_preprocessing(cfg);
    % calculate average
    cfg = [];
    cfg.channel = {'MEG'};
    avg3 = ft_timelockanalysis(cfg,data3);
    
    % plot the results for comparison
    tMEG = avg1.time;
    figure
    plot(tMEG, avg1.avg(1,:),'b')
    set(gca,'fontsize',14);
    hold on
    plot(tMEG, avg2.avg(1,:),'g')
    plot(tMEG, avg3.avg(1,:),'r')
    legend('avg1', 'avg2','avg3');
    plot(tMEG, avg1.avg,'b')
    hold on
    plot(tMEG, avg2.avg,'g')
    plot(tMEG, avg3.avg,'r')
    
    title('compare HB clean all channels')
    xlabel('Time (sec.)')
    ylabel('Amplitude (Tesla)');
    grid on
    
    if saveAvgs
        
        save avgs avg1 avg2 avg3
        
    end
    if saveFigs
        saveas(gcf, 'butterAll.fig');
        saveas(gcf, 'butterAll.png');
    end
end



