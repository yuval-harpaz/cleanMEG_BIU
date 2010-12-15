%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read Process Step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function step = read_proc_step(fid)

%all structures always start at byte sizeof(double)*N,
%where N is integer and sizeof(double) == 8
%(see <libdftk>/dftk_misc.C: int dftk_align(FILE *fp))
align_file(fid);

step.hdr.nbytes = fread(fid, 1, 'uint32=>uint32');
         type = char(fread(fid, 20, 'uchar'))';
%remove zeros from string
step.hdr.type = type(type>0);
step.hdr.checksum = fread(fid, 1, 'int32=>int32');

switch step.hdr.type
    case {'b_filt_hp' ...
          'b_filt_lp' ...
          'b_filt_notch'}
        step.frequency = fread(fid, 1, 'float32=>float32');
        step.reserved = fread(fid, 32, 'uchar')';
    case {'b_filt_b_pass' ...
          'b_filt_b_reject'}
        step.high_frequency = fread(fid, 1, 'float32=>float32');
        step.low_frequency = fread(fid, 1, 'float32=>float32');
        step.reserved = fread(fid, 32, 'uchar')';
        fseek(fid, 4, 'cof');%alignment
    otherwise   %user process 
        step.user_space_size = fread(fid, 1, 'uint32=>uint32');
        step.reserved = fread(fid, 32, 'uchar')';
        %read user data as array of bytes
        step.user_data = fread(fid, double(step.user_space_size), 'uchar');
        fseek(fid, mod(double(step.user_space_size),8), 'cof');%alignment
end