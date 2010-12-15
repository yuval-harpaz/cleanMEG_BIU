%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Process (dftk_proc_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function proc = read_proc_data(fid)

%all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) == 8
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

%alignment checked
proc.hdr.nbytes = fread(fid, 1, 'uint32=>uint32');
         type = char(fread(fid, 20, 'uchar'))';
proc.hdr.type = type(type>0);
proc.hdr.checksum = fread(fid, 1, 'int32=>int32');
     user = char(fread(fid, 32, 'uchar'))';
proc.user = user(user>0);
proc.timestamp = fread(fid, 1, 'uint32=>uint32');
     filename = char(fread(fid, 256, 'uchar'))';
proc.filename = filename(filename>0);
proc.total_steps = fread(fid, 1, 'uint32=>uint32');
proc.reserved = fread(fid, 32, 'uchar')';
fseek(fid, 4, 'cof');%alignment

%read process steps
for proc_step = 1:proc.total_steps
    proc.step{proc_step} = read_proc_step(fid);
end