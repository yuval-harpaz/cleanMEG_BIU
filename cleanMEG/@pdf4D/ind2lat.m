function [latency, epoch] = ind2lat(obj, index)

% IND2LAT [latency, epoch] = ind2lat(obj, index)
% converts data index into epoch number and latency

header = get(obj, 'header');

%total number of time point in pdf
total_points = pts_in_pdf(header);
index(index < 1) = 1;
index(index > total_points) = total_points;

if nargin < 2
    error('Too few arguments');
end

epoch = zeros(size(index));

for ep = 1:header.header_data.total_epochs
    %number of time point in the epoch
    lat2skip = double(header.epoch_data{ep}.pts_in_epoch);
    this_epoch = index <= lat2skip & epoch == 0;
    other_epoch = index > lat2skip & epoch == 0;
    epoch(this_epoch) = ep;
    index(other_epoch) = index(other_epoch) - lat2skip;
end

latency = double(index * header.header_data.sample_period - ...
    trigger_latency(header));