function [trig_start, trig_end] = trigger_latency(header)

%return trigger start and stop latency

trig_start = [];
trig_end = [];

for ev = 1:header.header_data.total_events
    if strcmp(header.event_data{ev}.event_name, 'Trigger')
        trig_start = header.event_data{ev}.start_lat;
        trig_end = header.event_data{ev}.end_lat;
        break
    end
end