function config = sysconf(varargin)

%Purpose:
% - read system config file

if nargin > 1
    filename = varargin{2};
    [p,f,ext] = fileparts(filename);
    if isempty(ext) || ~strcmp(ext, '.config')
        filename = sprintf('%s.config', filename);
    end
    if filename(1)~=filesep %this is relative filename
        %get STAGE
        [stat, stage] = unix('echo -n $STAGE');

        %no STAGE -> return
        if stat~=0
            config = [];
            return
        end
        filename = fullfile('/home', stage, 'config', filename);
    end
else
    %may be open listdlg to prompt for config name ?
    error('Too few inputs')
end

config = read_config(filename);