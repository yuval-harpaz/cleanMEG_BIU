function dataN = get_data_drives

%if something is wrong return empty cell array
dataN = {};

%get data path ('/home/whsbti/data')
data_path = get_data_path;
if isempty(data_path)
    return
end

%get data drives ('comeg1_data0' and so on)
tmp_data = dir(fullfile(data_path, '*data*'));
d = 0; %data drive index
for di=1:numel(tmp_data)
    if ~tmp_data(di).isdir
        continue
    end
    d = d +1;
    dataN{d} = fullfile(data_path, tmp_data(di).name);
end