%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read User Block Data (dftk_user_block_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function user_block_data = read_user_blok_data(fid)

%header and all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) is from C code
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

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