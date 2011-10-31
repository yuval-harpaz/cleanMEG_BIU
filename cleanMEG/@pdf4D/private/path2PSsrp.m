function id = path2PSsrp(pdf_path)

[run_path, pdf, ext] = fileparts(pdf_path);
id{5} = [pdf, ext];
[session_path, run, ext] = fileparts(run_path);
id{4} = [run, ext];
[scan_path, session, ext] = fileparts(session_path);
id{3} = dir2session([session, ext]);
[patient_path, scan, ext] = fileparts(scan_path);
id{2} = [scan, ext];
[data_path, patient, ext] = fileparts(patient_path);
id{1} = [patient, ext];