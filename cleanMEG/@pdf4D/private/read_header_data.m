%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read PDF Header (dftk_header_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = read_header_data(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is 8
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

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