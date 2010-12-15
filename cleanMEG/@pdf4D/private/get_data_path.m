function data_path = get_data_path

%get STAGE
[stat, stage] = unix('echo -n $STAGE');

%no STAGE -> return
if stat~=0
    data_path = '';
    return
end

data_path = fullfile('/home', stage, 'data');
