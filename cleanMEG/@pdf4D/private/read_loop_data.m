%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Meg Loop Data (dftk_loop_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loop_data = read_loop_data(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is from C code
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

loop_data.position = fread(fid, 3, 'double');
loop_data.direction = fread(fid, 3, 'double');
loop_data.radius = fread(fid, 1, 'double');
loop_data.wire_radius = fread(fid, 1, 'double');
loop_data.turns = fread(fid, 1, 'uint16=>uint16');
    fseek(fid, 2, 'cof');%alignment
loop_data.checksum = fread(fid, 1, 'int32=>int32');
loop_data.reserved = fread(fid, 32, 'uchar=>uchar')';