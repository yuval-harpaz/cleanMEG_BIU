function pdf_path = post2path(post)

%convert post string into pdf_path

%if something is wrong return empty string
pdf_path = '';

%get data drives ('/home/whsbti/data/comeg1_data0' and so on)
dataN = get_data_drives;
if isempty(dataN)
    return
end

%convert long string into cell array
%({patien_id, scan_id, session_id, run_id, pdf_id})
id = post2PSsrp(post);

%for pdf numel(id) == 5 (patient, sacn, session, run, pdf)
if length(id)<5
    return
end

%convert session id into directory name
id{3} = session2dir(id{3});

%search data drives for posted pdf
for di=1:numel(dataN)
    full_path = fullfile(dataN{di}, ...
        id{1}, id{2}, id{3}, id{4}, id{5});
    if exist(full_path, 'file')==2
        pdf_path = full_path;
    end
end