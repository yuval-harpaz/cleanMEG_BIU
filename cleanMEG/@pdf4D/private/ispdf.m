function test = ispdf(filename)

%test = ispdf(filename)
%
%Test MSI pdf header
%
%Revision 1.0  09/11/06  eugene.kronberg@uchsc.edu

if nargin ~= 1
    error('Wrong number of input arguments');
end

if ~exist(filename, 'file')
    test = false;
    return
end

%always big endian (Sun or Linux)
fid = fopen(filename, 'r', 'b');

if fid == -1
    error('Cannot open file %s', filename);
end

%last 8 bytes of the pdf is header offset
fseek(fid, -8, 'eof');
header_offset = fread(fid,1,'uint64');

%test for header offset could be here

%first byte of the header
fseek(fid, header_offset, 'bof');

%read dftk_header_data
header.header_data = read_header_data(fid);

fclose(fid);

%check if it's real pdf, otherwise return false
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

        test = false;
        return
end

%looks as filename IS pdf
test = true;