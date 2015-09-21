% this is how we deal with subjects who will not allow their hearbeat to be
% cleaned in the normal way

% This is what we do:
% 1)preprocess the data as one long trl in order to run ICA and get the HB trace (ECG)
% 2) Then insert it into the yuval' HBcorrect function which cleans the data.

% What we need: our raw data or the product of Idan's cleaning (minus HB cleaning)

% preprocess data into one long trl
source='xc,lf_c,rfhp0.1Hz';
cfg=[];
cfg.dataset=source;
cfg.trialfun='trialfun_raw';
cfg1=ft_definetrial(cfg);
cfg1.channel='MEG'; 
cfg1.bpfilter = 'yes';
cfg1.bpfreq = [1 30];
dataorig=ft_preprocessing(cfg1);

save dataorig dataorig -v7.3 %(this addition is necessary for saving very long files)

% down sample data
cfg            = [];
cfg.resamplefs = 100;
cfg.detrend    = 'no';
downSamp       = ft_resampledata(cfg, dataorig);

% run ica
cfg = [];
cfg.channel = {'MEG'};
findECG = ft_componentanalysis(cfg, downSamp);

save findECG findECG 

%%
% clear; close all;
% load findECG; load dataorig;

% see the components and find the HB component
cfg=[];
cfg.layout='4D248.lay';
cfg.channel = {findECG.label{1:10}};
cfg.continuous='yes';
cfg.viewmode='component';
cfg.blocksize=2;
comppic=ft_databrowser(cfg,findECG);

HBcomp = ?; % insert manually the component you think is the HB

%run the ICA in the original data
topo=findECG.topo;
topolabel=findECG.topolabel;
cfg = [];
cfg.topo      = topo;
cfg.topolabel = topolabel;
comp          = ft_componentanalysis(cfg, dataorig);
HBtrace = comp.trial{1,1}(HBcomp,:);
% figure; plot(HBtrace); % Here we can check the quality of the HB trace

save HBtrace HBtrace  

%% re-cleaning the data with the HBtrace
   
newData=correctHB([],[],[],HBtrace,[]); 
rewrite_pdf(newData,[],[],'xc,hb,lf'); 










