function chan_type = channel_type(obj, chanindx)

% CHANNEL_TYPE chan_type = channel_name(obj, chanindx)
% returns cell array of channel types
% chanindx - indecies of the channels

%BTi channel types
%(lower case, so that 'trigger' is type and 'TRIGGER' is channel)
bti_type = lower({ ...
    'MEG' ...
    'EEG' ...
    'REFERENCE' ...
    'EXTERNAL' ...
    'TRIGGER' ...
    'UTILITY' ...
    'DERIVED' ...
    'SHORTED' ...
    });

header = get(obj, 'header');
config = get(obj, 'config');
for ch = 1:length(chanindx)
    chan_no = header.channel_data{chanindx(ch)}.chan_no;
    chan_type{ch} = bti_type(config.channel_data{chan_no}.type);
end