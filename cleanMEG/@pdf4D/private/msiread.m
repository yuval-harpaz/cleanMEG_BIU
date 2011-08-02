function config = read_config(filename)

%config = read_config(filename)
%
%Reads MSI config file
%
%Converts C-style strings into MATLAB string by removing zeros
%Otherwise there is no data conversion

%Revision 1.0  09/11/06  eugene.kronberg@uchsc.edu

if nargin < 1
    error('Too few input argumets');
end

if ~exist(filename, 'file')
    error('File %s does not exist', filename);
end

fid = fopen(filename, 'r', 'b');

if fid == -1
    error('Cannot open file %s', filename);
end

%if we are not going to read data, or file is not pdf
data = [];

%defaul file type
if nargin < 2 | isempty(filetype)
    filetype = 'pdf';
end

switch filetype
    case 'pdf'
        if nargin < 3
            options = [];
        end
        if nargout > 1
            [header, data] = readpdf(fid, options);
        else
            header = readpdf(fid, options);
        end
    case 'config'
        if nargout > 1
            [header, data] = readconfig(fid);
        else
            header = readconfig(fid);
        end
    case 'hs_file'
        if nargout > 1
            [header, data] = readhsfile(fid);
        else
            header = readhsfile(fid);
        end
    otherwise
        fclose(fid);
        error('Wrong file type %s', filetype);
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Head Shape File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header, hs_points] = readhsfile(fid)

header.hdr.version = fread(fid, 1, 'uint32=>uint32');
header.hdr.timestamp = fread(fid, 1, 'int32=>int32');
header.hdr.checksum = fread(fid, 1, 'int32=>int32');
header.hdr.npoints = fread(fid, 1, 'int32=>int32');
header.lpa = fread(fid, 3, 'double');
header.rpa = fread(fid, 3, 'double');
header.nasion = fread(fid, 3, 'double');
header.cz = fread(fid, 3, 'double');
header.inion = fread(fid, 3, 'double');

if nargout > 1
    hs_points = fread(fid, [3 double(header.hdr.npoints)], 'double');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Config
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [config, ch_pos_n_dir] = readconfig(fid)

config.config_data = readconfigdata(fid);
config.Xfm = fread(fid, [4 4], 'double');
for ub = 1:config.config_data.total_user_blocks
    config.user_block_data{ub} = readuserblokdata(fid);
    %skip user block for now
    user_space_size = double(config.user_block_data{ub}.user_space_size);
    modulus = mod(user_space_size,8);
    if modulus
        padding = 8 - modulus;
    else
        padding = 0;
    end
    offset = user_space_size + padding;
    fseek(fid, offset, 'cof');
end
for ch = 1:config.config_data.total_chans
    config.channel_data{ch} = readchanneldata(fid);
end

if nargout > 1
    %positions and directions for the first meg loops
    ch_pos_n_dir = getchannelpos(config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Get Channel Position
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ch_pos_n_dir = getchannelpos(config)

meg_indx = [];
meg_ch = 0;
total_loops = [];

for ch = 1:config.config_data.total_chans
    name = config.channel_data{ch}.name;
    if config.channel_data{ch}.type == 1
        meg_ch = meg_ch + 1;
        meg_indx = [meg_indx str2num(name(2:end))];
        if isempty(total_loops)
            total_loops = config.channel_data{ch} ...
                .device_data.total_loops;%must be the same for all meg
        end
        for loop=1:total_loops
            ch_pos(meg_ch,1:3,loop) = ...
            config.channel_data{ch}.device_data.loop_data{loop}.position;
            ch_dir(meg_ch,1:3,loop) =  ...
            config.channel_data{ch}.device_data.loop_data{loop}.direction;
        end
        %config.channel_data{ch}.device_data
        [dummy indx] = sort(meg_indx);
    end
end

ch_pos = ch_pos(indx,:,:);
ch_dir = ch_dir(indx,:,:);

ch_pos_n_dir.position = ch_pos;
ch_pos_n_dir.direction = ch_dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Config Data (dftk_config_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function config_data = readconfigdata(fid)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read User Block Data (dftk_user_block_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function user_block_data = readuserblokdata(fid)

user_block_data.hdr.nbytes = fread(fid, 1, 'uint32=>uint32');
                    type = char(fread(fid, 20, 'uchar'))';
user_block_data.hdr.type = type(type>0);
user_block_data.hdr.checksum = fread(fid, 1, 'int32=>int32');
                user = char(fread(fid, 32, 'uchar'))';
user_block_data.user = user(user>0);
user_block_data.timestamp = fread(fid, 1, 'uint32=>uint32');
user_block_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
user_block_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Channel Data (dftk_channel_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channel_data = readchanneldata(fid)

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
        channel_data.device_data = readmegdevice(fid);
    case 2%eeg
        channel_data.device_data = readeegdevice(fid);
    case 4%external
        channel_data.device_data = readexternaldevice(fid);
    case 5%TRIGGER
        channel_data.device_data = readtriggerdevice(fid);
    case 6%utility
        channel_data.device_data = readutilitydevice(fid);
    case 7%derived
        channel_data.device_data = readderiveddevice(fid);
    case 8%shorted
        channel_data.device_data = readshorteddevice(fid);
    otherwise
        error('Unknown device type: %d\n', channel_data.type);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Meg Device Data (dftk_meg_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readmegdevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.inductance = fread(fid, 1, 'float32=>float32');
    fseek(fid, 4, 'cof');%alignment
device_data.Xfm = fread(fid, [4 4], 'double');
device_data.xform_flag = fread(fid, 1, 'uint16=>uint16');
device_data.total_loops = fread(fid, 1, 'uint16=>uint16');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

for loop = 1:device_data.total_loops
    device_data.loop_data{loop} = readloopdata(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Meg Loop Data (dftk_loop_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loop_data = readloopdata(fid)

loop_data.position = fread(fid, 3, 'double');
loop_data.direction = fread(fid, 3, 'double');
loop_data.radius = fread(fid, 1, 'double');
loop_data.wire_radius = fread(fid, 1, 'double');
loop_data.turns = fread(fid, 1, 'uint16=>uint16');
    fseek(fid, 2, 'cof');%alignment
loop_data.checksum = fread(fid, 1, 'int32=>int32');
loop_data.reserved = fread(fid, 32, 'uchar=>uchar')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Eeg Device Data (dftk_eeg_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readeegdevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.impedance = fread(fid, 1, 'float32=>float32');
    fseek(fid, 4, 'cof');%alignment
device_data.Xfm = fread(fid, [4 4], 'double');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read External Device Data (dftk_external_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readexternaldevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Trigger Device Data (dftk_trigger_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readtriggerdevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Utility Device Data (dftk_utility_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readutilitydevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Utility Device Data (dftk_derived_device_data ?)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readderiveddevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Shorted Device Data (dftk_shorted_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = readshorteddevice(fid)

device_data.hdr = readdeviceheader(fid);
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Device Header (dftk_device_header)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = readdeviceheader(fid)

header.size = fread(fid, 1, 'uint32=>uint32');
header.checksum = fread(fid, 1, 'int32=>int32');
header.reserved = fread(fid, 32, 'uchar=>uchar')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read PDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header, data] = readpdf(fid, options)

header = [];
data = [];

%last 8 bytes of the pdf is header offset
fseek(fid, -8, 'eof');
header_offset = fread(fid,1,'uint64');

%first byte of the header
fseek(fid, header_offset, 'bof');
header = readpdfhdr(fid);

%first byte of the data
fseek(fid, 0, 'bof');
if nargout > 1
    [header, data] = readpdfdata(fid, header, options);
else
    header = readpdfdata(fid, header, options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read PDF Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header, data] = readpdfdata(fid, header, options)

data = [];

%BTi Data Formats:
SHORT   =	1;
LONG    =	2;
FLOAT   =	3;
DOUBLE  =	4;

if isempty(options)
    return
end

%total number of channels in pdf (get before sort_channels)
total_chans = double(header.header_data.total_chans);

[header, chan_index] = sort_channels(header, options);

[header, lat_index] = sort_latency(header, options);

%no channels or no latencies to read
if nargout < 2 || isempty(chan_index) || isempty(lat_index)
    return
end

%set order of dimensions
if isfield(options, 'first_dim') && strcmp(options.first_dim, 'channel')
    chan_first = true;
    data = zeros(numel(chan_index), sum(lat_index));
else
    chan_first = false;
    data = zeros(sum(lat_index), numel(chan_index));
end

switch header.header_data.data_format
    case SHORT
        data_format = 'int16=>int16';
        time_slice_size = 2 * total_chans;
        data = int16(data);
    case LONG
        data_format = 'int32=>int32';
        time_slice_size = 4 * total_chans;
        data = int32(data);
    case FLOAT
        data_format = 'float32=>float32';
        time_slice_size = 4 * total_chans;
        data = single(data);
    case DOUBLE
        data_format = 'double';
        time_slice_size = 8 * total_chans;
        data = double(data);
    otherwise
        error('Wrong data format : %d\n', header.header_data.data_format);
end

lat_out = 0;
for lat = 1:numel(lat_index)
    if lat_index(lat)
        %read time slice
        time_slice = fread(fid, double(total_chans), data_format);
        lat_out = lat_out + 1;
        if chan_first
            data(:,lat_out) = time_slice(chan_index)';
        else
            data(lat_out,:) = time_slice(chan_index);
        end
    else
        %skip time slice
        fseek(fid, time_slice_size, 'cof');
    end
end

if chan_first
    %data dimensions: channel, time, epoch
    new_dims = [size(data,1), ...
    size(data,2) / header.header_data.total_epochs, ...
    header.header_data.total_epochs];
else
    %data dimensions: time, epoch, channel
    new_dims = [size(data,1) / header.header_data.total_epochs, ...
        header.header_data.total_epochs, size(data,2)];
end
data = reshape(data,new_dims);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sort Latency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header, lat_index] = sort_latency(header, options)

%get first latency from event 'Trigger'
trigger = trigger_start(header);

%count total number of time-points
%get first time-point for all epochs
epoch_start = ones(1,header.header_data.total_epochs+1);
for epoch = 1:header.header_data.total_epochs
    epoch_start(epoch+1) = epoch_start(epoch) + ...
        header.epoch_data{epoch}.pts_in_epoch;
end
epoch_end = epoch_start(2:end) - 1;
total_points = epoch_end(end);
epoch_start(end) = [];

lat_index = zeros(1, total_points);

if ~isfield(options, 'latency_list')
    %default all epochs, and all time-points
    lat_index(:) = 1;
    return
elseif isstruct(options.latency_list)
    latency_list{1} = options.latency_list;
elseif ~iscell(options.latency_list) | ...
       ~isfield(options.latency_list{1}, 'epochs') | ...
       ~isfield(options.latency_list{1}, 'min_lat') | ...
       ~isfield(options.latency_list{1}, 'max_lat')
       error('Wrong latency list\n');
else
    latency_list = options.latency_list;
end

for lat_grp = 1:numel(latency_list)
    if ischar(latency_list{lat_grp}.epochs)
        if strcmp(latency_list{lat_grp}.epochs, 'all')
            epochs = 1:header.header_data.total_epochs;
        else
            error('Wrong epoch list %s\n', ...
                latency_list{lat_grp}.epochs);
        end
    else
        epochs = latency_list{lat_grp}.epochs;
    end
    epochs(epochs < 1) = [];
    epochs(epochs > header.header_data.total_epochs) = [];
    if isempty(epochs)
        continue
    end
    for epoch = epochs
        ep_lat = ...
            double([1:header.epoch_data{epoch}.pts_in_epoch]) * ...
            header.header_data.sample_period - trigger;
        lat_index(epoch_start(epoch):epoch_end(epoch)) = ...
            lat_index(epoch_start(epoch):epoch_end(epoch)) | ...
            (ep_lat > latency_list{lat_grp}.min_lat & ...
             ep_lat < latency_list{lat_grp}.max_lat);
    end
end

for epoch = double(header.header_data.total_epochs):-1:1
    pts_in_epoch = ...
    uint32(sum(lat_index(epoch_start(epoch):epoch_end(epoch))));
    if pts_in_epoch
        header.epoch_data{epoch}.pts_in_epoch = pts_in_epoch;
        header.epoch_data{epoch}.epoch_duration = ...
            single(double(pts_in_epoch -1) * ...
                   double(header.header_data.sample_period));
    else
        %no points - delete epoch from the header
        header.epoch_data(epoch) = [];
    end
end

%set new total number of epochs
header.header_data.total_epochs = uint32(numel(header.epoch_data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Trigger Start
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trigger = trigger_start(header)

%if there is no trigger (for whatever reason) start from zero
trigger = 0;

if ~isfield(header, 'event_data')
    return
end

%find 'Trigger' event
for event =1:numel(header.event_data)
    if strcmp(header.event_data{event}.event_name, 'Trigger')
        trigger = header.event_data{event}.start_lat;
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sort Channel List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header, chan_index] = sort_channels(header, options)

chan_index = [];

if ~isfield(options, 'channel_list')
    chan_index = [1:header.header_data.total_chans];
    return
elseif isempty(options.channel_list)
    error('Empty channel list\n');
elseif ischar(options.channel_list)
    switch options.channel_list
        case 'all'
            chan_index = [1:header.header_data.total_chans];
            return
        case {'all_meg' 'meg'}
            meg_ch = 0;
            chan_index = [];
            %skip channel which start with A but not meg (like 'Accel X')
            %meg channel should have number after first A
            for ch = 1:header.header_data.total_chans
                if ~(header.channel_data{ch}.chan_label(1)=='A') | ...
                   isempty(str2num(header.channel_data{ch}.chan_label(2:end)))
                    continue
                end
                meg_ch = meg_ch + 1;
                chan_index = [chan_index ch];
                meg_indx(meg_ch) = ...
                    str2num(header.channel_data{ch}.chan_label(2:end));
                meg_lbl{meg_ch} = header.channel_data{ch}.chan_label;
                %meg_no(meg_ch) = header.channel_data{ch}.chan_no;
            end
            [dummy, index] = sort(meg_indx);
            chan_index = chan_index(index);
            header.channel_data = header.channel_data(chan_index);
            return
        case 'left_meg'
            error('Not implemented yet :( Channel list: %s\n', ...
                options.channel_list);
        case 'right_meg'
            error('Not implemented yet :( Channel list: %s\n', ...
                options.channel_list);
        case 'eeg'
            eeg_ch = 0;
            for ch = 1:header.header_data.total_chans
                if ~(header.channel_data{ch}.chan_label(1)=='E')
                    continue
                end
                eeg_ch = eeg_ch + 1;
                eeg_indx(eeg_ch) = ...
                    str2num(header.channel_data{ch}.chan_label(2:end));
                eeg_lbl{eeg_ch} = header.channel_data{ch}.chan_label;
                eeg_no(eeg_ch) = header.channel_data{ch}.chan_no;
            end
            [dummy, index] = sort(eeg_indx);
            chan_index = eeg_no(index);
            header.channel_data = header.channel_data(chan_index);
            return
        case 'ext'
            ext_ch = 0;
            for ch = 1:header.header_data.total_chans
                if ~(header.channel_data{ch}.chan_label(1)=='X')
                    continue
                end
                ext_ch = ext_ch + 1;
                ext_indx(ext_ch) = ...
                    str2num(header.channel_data{ch}.chan_label(2:end));
                ext_lbl{ext_ch} = header.channel_data{ch}.chan_label;
                ext_no(ext_ch) = header.channel_data{ch}.chan_no;
            end
            [dummy, index] = sort(ext_indx);
            chan_index = ext_no(index);
            header.channel_data = header.channel_data(chan_index);
            return
        case {'all_ref' 'ref'}
            ref_ch = 0;
            for ch = 1:header.header_data.total_chans
                if ~(header.channel_data{ch}.chan_label(1)=='M' | ...
                     header.channel_data{ch}.chan_label(1)=='G')
                    continue
                end
                ref_ch = ref_ch + 1;
                chan_index = [chan_index ch];
            end
            header.channel_data = header.channel_data(chan_index);
            return
        case 'g_ref'
            ref_ch = 0;
            for ch = 1:header.header_data.total_chans
                if ~(header.channel_data{ch}.chan_label(1)=='G')
                    continue
                end
                ref_ch = ref_ch + 1;
                chan_index = [chan_index ch];
            end
            header.channel_data = header.channel_data(chan_index);
            return
        case 'm_ref'
            ref_ch = 0;
            for ch = 1:header.header_data.total_chans
                if ~(header.channel_data{ch}.chan_label(1)=='M')
                    continue
                end
                ref_ch = ref_ch + 1;
                chan_index = [chan_index ch];
            end
            header.channel_data = header.channel_data(chan_index);
            return
        otherwise
            error('Wrong channel list: %s\n', options.channel_list);
    end
end

if iscell(options.channel_list)
    chan_index = zeros(1,numel(options.channel_list));
    for hdr_ch = 1:header.header_data.total_chans
        chan_index(find(strcmp(header.channel_data{hdr_ch}.chan_label, ...
                    options.channel_list))) = hdr_ch;
    end
    chan_index(~chan_index) = [];
    header.channel_data = header.channel_data(chan_index);
else
    error('Wrong channel list field\n');
end

%set new total_chans
header.header_data.total_chans = uint16(numel(chan_index));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read PDF Header
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = readpdfhdr(fid)

%use cell arrays (for epochs and procs) since structures could be different
%use cell arrays for the rest, just to be consistent

%read dftk_header_data
header.header_data = readheader(fid);

%check if it's real pdf, otherwise return []
if isempty(header.header_data.version) || ...
   isempty(header.header_data.file_type) || ...
   isempty(header.header_data.data_format) || ...
   isempty(header.header_data.acq_mode) || ...
   isempty(header.header_data.total_epochs) || ...
   isempty(header.header_data.input_epochs) || ...
   isempty(header.header_data.total_events) || ...
   isempty(header.header_data.total_fixed_events) || ...
   isempty(header.header_data.sample_period) || ...
   isempty(header.header_data.xaxis_label) || ...
   isempty(header.header_data.total_processes) || ...
   isempty(header.header_data.total_chans) || ...
   isempty(header.header_data.checksum) || ...
   isempty(header.header_data.total_ed_classes) || ...
   isempty(header.header_data.total_associated_files) || ...
   isempty(header.header_data.last_file_index) || ...
   isempty(header.header_data.timestamp) || ...
   isempty(header.header_data.reserved)

        header = [];
        return
end

%read dftk_epoch_data
for epoch = 1:header.header_data.total_epochs;
    header.epoch_data{epoch} = readepoch(fid);
end

%read dftk_channel_ref_data
for channel = 1:header.header_data.total_chans
    header.channel_data{channel} = readchannelref(fid);
end

%read dftk_event_data
for event = 1:header.header_data.total_fixed_events
    header.event_data{event} = readevent(fid);
end

%it might not work for all processes, we stop here
%(to try just comment next line)
return

%read dftk_proc_data
for proc = 1:header.header_data.total_processes
    header.proc_data{proc} = readprocess(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Header (dftk_header_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = readheader(fid)

%alignment checked
header.version = fread(fid, 1, 'uint16=>uint16');
       file_type = char(fread(fid, 5, 'uchar'))';
header.file_type = file_type(file_type>0);
fseek(fid, 1, 'cof');%alignment
header.data_format = fread(fid, 1, 'int16=>int16');
header.acq_mode = fread(fid, 1, 'uint16=>uint16');
header.total_epochs = fread(fid, 1, 'uint32=>uint32');
header.input_epochs = fread(fid, 1, 'uint32=>uint32');
header.total_events = fread(fid, 1, 'uint32=>uint32');
header.total_fixed_events = fread(fid, 1, 'uint32=>uint32');
header.sample_period = fread(fid, 1, 'float32=>float32');
       xaxis_label = char(fread(fid, 16, 'uchar'))';
header.xaxis_label = xaxis_label(xaxis_label>0);
header.total_processes = fread(fid, 1, 'uint32=>uint32');
header.total_chans = fread(fid, 1, 'uint16=>uint16');
fseek(fid, 2, 'cof');%alignment
header.checksum = fread(fid, 1, 'int32=>int32');
header.total_ed_classes = fread(fid, 1, 'uint32=>uint32');
header.total_associated_files = fread(fid, 1, 'uint16=>uint16');
header.last_file_index = fread(fid, 1, 'uint16=>uint16');
header.timestamp = fread(fid, 1, 'uint32=>uint32');
header.reserved = fread(fid, 20, 'uchar')';
fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Epoch (dftk_epoch_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epoch = readepoch(fid)

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
    epoch.var_event{event} = readevent(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Channel Reference (dftk_channel_ref_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channel = readchannelref(fid)

%alignment checked
        chan_label = char(fread(fid, 16, 'uchar'))';
channel.chan_label = chan_label(chan_label>0);
channel.chan_no = fread(fid, 1, 'uint16=>uint16');
channel.attributes = fread(fid, 1, 'uint16=>uint16');
channel.scale = fread(fid, 1, 'float32=>float32');
        yaxis_label = char(fread(fid, 16, 'uchar'))';
channel.yaxis_label = yaxis_label(yaxis_label>0);
channel.valid_min_max = fread(fid, 1, 'uint16=>uint16');
fseek(fid, 6, 'cof');%alignment
channel.ymin = fread(fid, 1, 'float64');
channel.ymax = fread(fid, 1, 'float64');
channel.index = fread(fid, 1, 'uint32=>uint32');
channel.checksum = fread(fid, 1, 'int32=>int32');
%something new?
channel.whatisit = char(fread(fid, 4, 'uchar'))';
channel.reserved = fread(fid, 28, 'uchar')';
%reserved first 4 bytes are not 0?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Event (dftk_event_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function event = readevent(fid)

%alignment checked
      event_name = char(fread(fid, 16, 'uchar'))';
event.event_name = event_name(event_name>0);
event.start_lat = fread(fid, 1, 'float32=>float32');
event.end_lat = fread(fid, 1, 'float32=>float32');
event.step_size = fread(fid, 1, 'float32=>float32');
event.fixed_event = fread(fid, 1, 'uint16=>uint16');
fseek(fid, 2, 'cof');%alignment
event.checksum = fread(fid, 1, 'int32=>int32');
event.reserved = fread(fid, 32, 'uchar')';
fseek(fid, 4, 'cof');%alignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Process (dftk_proc_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function proc = readprocess(fid)

%alignment checked
proc.hdr.nbytes = fread(fid, 1, 'uint32=>uint32');
         type = char(fread(fid, 20, 'uchar'))';
proc.hdr.type = type(type>0);
proc.hdr.checksum = fread(fid, 1, 'int32=>int32');
     user = char(fread(fid, 32, 'uchar'))';
proc.user = user(user>0);
proc.timestamp = fread(fid, 1, 'uint32=>uint32');
     filename = char(fread(fid, 256, 'uchar'))';
proc.filename = filename(filename>0);
proc.total_steps = fread(fid, 1, 'uint32=>uint32');
proc.reserved = fread(fid, 32, 'uchar')';
fseek(fid, 4, 'cof');%alignment

%read process steps
for proc_step = 1:proc.total_steps
    proc.step{proc_step} = readprocstep(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Process Step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function step = readprocstep(fid)

step.hdr.nbytes = fread(fid, 1, 'uint32=>uint32');
         type = char(fread(fid, 20, 'uchar'))';
%remove zeros from string
step.hdr.type = type(type>0);
step.hdr.checksum = fread(fid, 1, 'int32=>int32');

switch step.hdr.type
    case {'b_filt_hp' ...
          'b_filt_lp' ...
          'b_filt_notch'}
        step.frequency = fread(fid, 1, 'float32=>float32');
        step.reserved = fread(fid, 32, 'uchar')';
    case {'b_filt_b_pass' ...
          'b_filt_b_reject'}
        step.high_frequency = fread(fid, 1, 'float32=>float32');
        step.low_frequency = fread(fid, 1, 'float32=>float32');
        step.reserved = fread(fid, 32, 'uchar')';
        fseek(fid, 4, 'cof');%alignment
    otherwise   %user process 
        step.user_space_size = fread(fid, 1, 'uint32=>uint32');
        step.reserved = fread(fid, 32, 'uchar')';
end
