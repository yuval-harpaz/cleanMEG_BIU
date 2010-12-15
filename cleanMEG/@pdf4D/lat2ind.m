function index = lat2ind(obj, epoch, latency)

% LAT2IND index = lat2ind(obj, epoch, latency)
% converts epoch number and latency into data indecies

header = get(obj, 'header');

if nargin < 3
    error('Too few arguments');
end

if isempty(epoch) || epoch < 1
    epoch = 1;
elseif epoch > header.header_data.total_epochs
    epoch = header.header_data.total_epochs;
else
    epoch = round(epoch);
end

lat2skip = 0;
for ep = 1:epoch-1
    lat2skip = lat2skip + double(header.epoch_data{ep}.pts_in_epoch);
end

index = double(round( ...
    (latency + trigger_latency(header)) / ...
    header.header_data.sample_period)) + 1;
index(index > header.epoch_data{epoch}.pts_in_epoch) = ...
    header.epoch_data{epoch}.pts_in_epoch;
index(index < 1) = 1;
index = index + lat2skip;