%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Device Header (dftk_device_header)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = read_device_header(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is from C code
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

header.size = fread(fid, 1, 'uint32=>uint32');
header.checksum = fread(fid, 1, 'int32=>int32');
header.reserved = fread(fid, 32, 'uchar=>uchar')';
