%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Config Data (dftk_config_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function config_data = read_config_data(fid)

config_data.version = fread(fid, 1, 'uint16=>uint16');
            site_name = char(fread(fid, 32, 'uchar'))';
config_data.site_name = site_name(site_name>0);
            dap_hostname = char(fread(fid, 16, 'uchar'))';
config_data.dap_hostname = dap_hostname(dap_hostname>0);
config_data.sys_type = fread(fid, 1, 'uint16=>uint16');
config_data.sys_options = fread(fid, 1, 'uint32=>uint32');
config_data.supply_freq = fread(fid, 1, 'uint16=>uint16');
config_data.total_chans = fread(fid, 1, 'uint16=>uint16');
config_data.system_fixed_gain = fread(fid, 1, 'float32=>float32');
config_data.volts_per_bit = fread(fid, 1, 'float32=>float32');
config_data.total_sensors = fread(fid, 1, 'uint16=>uint16');
config_data.total_user_blocks = fread(fid, 1, 'uint16=>uint16');
config_data.next_derived_channel_number = fread(fid, 1, 'uint16=>uint16');
    fseek(fid, 2, 'cof');%alignment
config_data.checksum = fread(fid, 1, 'int32=>int32');
config_data.reserved = fread(fid, 32, 'uchar=>uchar')';