function session_dir = session2dir(session_id)

%convert session_id into directory name

session_dir = session_id;
session_dir(session_id=='/') = '%';
session_dir(session_id==' ') = '@';
