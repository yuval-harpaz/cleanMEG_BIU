function header = read_header(filename)

%header = read_header(filename)
%
%Reads MSI pdf header
%
%Converts C-style strings into MATLAB string by removing zeros
%Otherwise there is no data conversion

%Revision 1.0  09/11/06  eugene.kronberg@uchsc.edu

%Revision 1.1  11/08/06, header starts at byte 8*N, where N is integer

if nargin ~= 1
    error('Wrong number of input arguments');
end

if ~ispdf(filename)
    error('File %s is not pdf', filename);
end

%always big endian (Sun or Linux)
fid = fopen(filename, 'r', 'b');

if fid == -1
    error('Cannot open file %s', filename);
end

%last 8 bytes of the pdf is header offset
fseek(fid, -8, 'eof');
header_offset = fread(fid,1,'uint64');

%first byte of the header
%(read_header_data will call align_file for 8*N alignment)
fseek(fid, header_offset, 'bof');

%use cell arrays (for epochs and procs) since structures could be different
%use cell arrays for the rest, just to be consistent

%read dftk_header_data
header.header_data = read_header_data(fid);

%read dftk_epoch_data
for epoch = 1:header.header_data.total_epochs;
    header.epoch_data{epoch} = read_epoch_data(fid);
end

%read dftk_channel_ref_data
for channel = 1:header.header_data.total_chans
    header.channel_data{channel} = read_channel_ref_data(fid);
end

%read dftk_event_data
for event = 1:header.header_data.total_fixed_events
    header.event_data{event} = read_event_data(fid);
end

%it might not work for all processes, so we stop here
% %read dftk_proc_data
% for proc = 1:header.header_data.total_processes
%     header.proc_data{proc} = read_proc_data(fid);
% end

fclose(fid);