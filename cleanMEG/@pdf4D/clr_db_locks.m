function clr_db_locks(obj)

%run clr_db_locks on patient

%if there is no patient id return
if isempty(get(obj, 'patient'));
    fprintf('No Patient ID\n');
    return
end

%run clr_db_locks
[s,w]=unix(sprintf('clr_db_locks -P %s', get(obj, 'patient')));

%if there was an error, pass it to user
if s~=0
    waitfor(errordlg(w, 'clr_db_locks error'));
end