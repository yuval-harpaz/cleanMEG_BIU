function config = read_config(filename)

%config = read_config(filename)
%
%Reads MSI (system or run) config file
%
%Converts C-style strings into MATLAB string by removing zeros
%Otherwise there is no data conversion

%Revision 1.0  09/11/06  eugene.kronberg@uchsc.edu

if nargin < 1
    error('Too few input argumets');
end

if ~isconfig(filename)
    error('File %s is not valid config file', filename);
end

fid = fopen(filename, 'r', 'b');

if fid == -1
    error('Cannot open file %s', filename);
end

config.config_data = read_config_data(fid);

%alignment is fine
config.Xfm = fread(fid, [4 4], 'double');

%user blocks
for ub = 1:config.config_data.total_user_blocks
    config.user_block_data{ub} = read_user_blok_data(fid);
    %skip user block for now
    user_space_size = double(config.user_block_data{ub}.user_space_size);
    fseek(fid, user_space_size, 'cof');
    %old alignment (before function align_file was written)
%     modulus = mod(user_space_size,8);
%     if modulus
%         padding = 8 - modulus;
%     else
%         padding = 0;
%     end
%     offset = user_space_size + padding;
%     fseek(fid, offset, 'cof');
end

%channels
for ch = 1:config.config_data.total_chans
    config.channel_data{ch} = read_channel_data(fid);
end

fclose(fid);