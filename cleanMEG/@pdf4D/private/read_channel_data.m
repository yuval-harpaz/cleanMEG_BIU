%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Channel Data (dftk_channel_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channel_data = read_channel_data(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is from C code
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

             name = char(fread(fid, 16, 'uchar'))';
channel_data.name = name(name>0);
channel_data.chan_no = fread(fid, 1, 'uint16=>uint16');
channel_data.type = fread(fid, 1, 'uint16=>uint16');
channel_data.sensor_no = fread(fid, 1, 'int16=>int16');
    fseek(fid, 2, 'cof');%alignment
channel_data.gain = fread(fid, 1, 'float32=>float32');
channel_data.units_per_bit = fread(fid, 1, 'float32=>float32');
             yaxis_label = char(fread(fid, 16, 'uchar'))';
channel_data.yaxis_label = yaxis_label(yaxis_label>0);
channel_data.aar_val = fread(fid, 1, 'double');
channel_data.checksum = fread(fid, 1, 'int32=>int32');
channel_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

switch channel_data.type
    case {1, 3}%meg/ref
        channel_data.device_data = read_meg_device_data(fid);
    case 2%eeg
        channel_data.device_data = read_eeg_device_data(fid);
    case 4%external
        channel_data.device_data = read_external_device_data(fid);
    case 5%TRIGGER
        channel_data.device_data = read_trigger_device_data(fid);
    case 6%utility
        channel_data.device_data = read_utility_device_data(fid);
    case 7%derived
        channel_data.device_data = read_derived_device_data(fid);
    case 8%shorted
        channel_data.device_data = read_shorted_device_data(fid);
    otherwise
        error('Unknown device type: %d\n', channel_data.type);
end
