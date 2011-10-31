%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read External Device Data (dftk_external_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = read_external_device_data(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is from C code
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

device_data.hdr = read_device_header(fid);
device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');%alignment
