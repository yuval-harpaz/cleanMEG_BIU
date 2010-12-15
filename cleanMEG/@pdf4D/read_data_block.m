function data_block = read_data_block(obj, lat, chanindx)

% READ_DATA_BLOCK data = read_data_block(obj, latindx, chanindx)
% chanindx must be vector of channel indecies, default read all channels
% latindx could have
%    no elements ([]) to read all time points
%    one element (for one time slice)
%    two elemnts for min and max latency indecies

% Revision 0.1.1  02/16/07 eugene
% remove variable "data"

% Revision 0.1.4  03/12/08 eugene
% fix "short-format" bug

%that was before setting defaults (all argins should be present):
% %not enough argins
% if nargin < 3
%     error('Too few input argumets');
% end

%pdf header
header = get(obj, 'Header');

%empty header means no pdf
if isempty(header)
    error('Need pdf to read data from')
end

%total number of channels in pdf
total_chans = double(header.header_data.total_chans);

%channel indecies (if empty -> read all channels)
if nargin < 3 || isempty(chanindx)
    %default is to read all channels
    chanindx = 1:total_chans;
end

%latency indexes (if empty -> read all time points)
if nargin < 2
    lat = [];
end

switch length(lat)
    case 0
        %default is to read all time points
        lat2skip = 0;
        lat2read = double(pts_in_pdf(header));
    case 1
        %read just one time slice
        lat2skip = lat - 1;
        lat2read = 1;
    case 2
        %read data block
        lat2skip = lat(1) - 1;
        lat2read = lat(2) - lat2skip;
    otherwise
        error('Wrong latency vector length: %d', length(lat))
end

%init data matrix
% data = zeros(length(chanindx), lat2read);

%BTi Data Formats:
SHORT   =	1;
LONG    =	2;
FLOAT   =	3;
DOUBLE  =	4;

switch header.header_data.data_format
    case SHORT
        data_format = 'int16=>int16';
        time_slice_size = 2 * total_chans;
%         data = int16(data);
        scale_data = true;
    case LONG
        data_format = 'int32=>int32';
        time_slice_size = 4 * total_chans;
%         data = int32(data);
        scale_data = true;
    case FLOAT
        data_format = 'float32=>float32';
        time_slice_size = 4 * total_chans;
%         data = single(data);
        scale_data = false;
    case DOUBLE
        data_format = 'double';
        time_slice_size = 8 * total_chans;
%         data = double(data);
        scale_data = false;
    otherwise
        error('Wrong data format : %d\n', header.header_data.data_format);
end

%%%%%%%%%%%%%%%% reading data

%open file, always big-endean
fid = fopen(obj.FileName, 'r', 'b');

if fid == -1
    error('Cannot open file %s', obj.FileName);
end

%skip some time slices
num2skip = time_slice_size * round(lat2skip);
if num2skip>0
    status=fseek(fid, num2skip, 'bof');
    if status~=0
        error('pdf4D:readFile', 'Didnot skip the beginning')
    end
elseif num2skip<0
    error('MATLAB:pdf4D:wrongLatencies','Cannot start from negative latency')
end
data_block = fread(fid, double([total_chans, round(lat2read)]), data_format);

%select (and sort) channels
% data(:,:) = data_block(chanindx,:);
data_block = data_block(chanindx,:);

%close data file
fclose(fid);

%%%%%%%%%%%%%%%% scale channels if needed
if scale_data
    config = get(obj, 'config');
    if isempty(config)
        warning('pdf4D:DataScale', 'No Config: Need to scale data!\n');
        return
    end
    scale = channel_scale(config, header, chanindx);
%     data = single(data) .* repmat(scale', 1, size(data,2));
    data_block = single(data_block) .* repmat(scale', 1, size(data_block,2));
end