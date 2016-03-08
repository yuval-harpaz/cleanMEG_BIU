%% OMEGA database
% relies on ctf2matlab:
% http://kurage.nimh.nih.gov/library/Meg/ctf2matlab.tgz
% to save a FieldTrip structure you have to have FieldTrip too


%% use default options, save a fieldtrip structure to 'correctHBdata.mat;
% cd into the *.ds folder
correctHB;

%% Rtopo
% get the cleaned data as a variable, the times of the heartbeat and the
% topography at the peak
[cleaned,HBtimes,~,~,~,Rtopo]=correctHB;

