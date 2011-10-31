function chan_name = channel_name(obj, chanindx)

% CHANNEL_NAME chan_name = channel_name(obj, chanindx)
% returns cell array of channel names
% chanindx - indecies of the channels

header = get(obj, 'header');
config = get(obj, 'config');
for ch = 1:length(chanindx)
    chan_no = header.channel_data{chanindx(ch)}.chan_no;
    chan_name{ch} = config.channel_data{chan_no}.name;
end