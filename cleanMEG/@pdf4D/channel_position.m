function chan_pos = channel_position(obj, chanindx)

% CHANNEL_POSITION chan_pos = channel_position(obj, chanindx)
% returns array of structures with fields "position" and "direction"
%    for meg channels both fields are 3 by N, where N is number of coils
%    for non-meg channels both fields are empty
% chanindx - indecies of the channels

header = get(obj, 'header');
config = get(obj, 'config');

for ch = 1:length(chanindx)
    chan_no = header.channel_data{chanindx(ch)}.chan_no;
    switch config.channel_data{chan_no}.type
        case {1 3} %meg or ref
            nl = config.channel_data{chan_no}.device_data.total_loops;
            position = [];
            direction = [];
            for l = 1:nl
                position(:,l) = ...
                    config.channel_data{chan_no} ...
                    .device_data.loop_data{l}.position;
                direction(:,l) = ...
                    config.channel_data{chan_no} ...
                    .device_data.loop_data{l}.direction;
            end
            chan_pos(ch) = struct( ...
                'position', position, ...
                'direction', direction );
        otherwise
            chan_pos(ch) = struct( ...
                'position', [], ...
                'direction', [] );
    end
end