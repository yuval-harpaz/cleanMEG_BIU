function chan_label = channel_label(obj, chanindx)

% CHANNEL_LABEL chan_label = channel_label(obj, chanindx)
% returns cell array of channel labels
% chanindx - indecies of the channels

header = get(obj, 'header');
for ch = 1:length(chanindx)
    chan_label{ch} = header.channel_data{chanindx(ch)}.chan_label;
end