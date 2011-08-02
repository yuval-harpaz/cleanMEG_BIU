%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Epoch (dftk_epoch_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epoch = read_epoch_data(fid)

%all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) == 8
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

%alignment checked
epoch.pts_in_epoch = fread(fid, 1, 'uint32=>uint32');
epoch.epoch_duration = fread(fid, 1, 'float32=>float32');
epoch.expected_iti = fread(fid, 1, 'float32=>float32');
epoch.actual_iti = fread(fid, 1, 'float32=>float32');
epoch.total_var_events = fread(fid, 1, 'uint32=>uint32');
epoch.checksum = fread(fid, 1, 'int32=>int32');
epoch.epoch_timestamp = fread(fid, 1, 'int32=>int32');
epoch.reserved = fread(fid, 28, 'uchar')';

%read dftk_event_data (var_events)
for event = 1:epoch.total_var_events
    epoch.var_event{event} = read_event_data(fid);
end