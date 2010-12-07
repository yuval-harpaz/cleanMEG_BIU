%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Channel Reference (dftk_channel_ref_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channel = read_channel_ref_data(fid)

%all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) == 8
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

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