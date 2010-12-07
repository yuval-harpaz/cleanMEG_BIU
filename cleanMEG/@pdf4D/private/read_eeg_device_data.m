%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Eeg Device Data (dftk_eeg_device_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function device_data = read_eeg_device_data(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is from C code
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

device_data.hdr = read_device_header(fid);
device_data.impedance = fread(fid, 1, 'float32=>float32');
    fseek(fid, 4, 'cof');%alignment
device_data.Xfm = fread(fid, [4 4], 'double');
device_data.reserved = fread(fid, 32, 'uchar=>uchar')';