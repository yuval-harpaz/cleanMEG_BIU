function scale = channel_scale(config, header, chanindx)

%channel scale (single, same as units_per_bit)
scale = single(zeros(1, length(chanindx)));

chan_no = zeros(1, length(chanindx));
%sort channel structures
for ch = 1:numel(chanindx)
    chan_no = header.channel_data{chanindx(ch)}.chan_no;
    scale(ch) = ...
        config.channel_data{chan_no}.units_per_bit;
%         config.channel_data{chan_no}.gain;
end