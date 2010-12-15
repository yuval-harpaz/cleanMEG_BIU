function session_id = dir2session(session_dir)

%convert session_id into directory name

session_id = session_dir;
session_id(session_dir=='%') = '/';
session_id(session_dir=='@') = ' ';
