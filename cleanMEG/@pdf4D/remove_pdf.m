function remove_pdf(obj, pdf_id)

% REMOVE_PDF remove_pdf(obj, pdf_id)
% pdf_id - filename(s) for the pdf (the rest of ids taken from the obj)
%          for multiple files should be cell array

if ischar(pdf_id)
    pdf_id = {pdf_id};
end

for ii=1:numel(pdf_id)
    %command to run asc_to_pdf
    run_remover = sprintf( ...
        'remover -P "%s" -S "%s" -s "%s" -r "%s" -p "%s"\n', ...
        get(obj, 'patient'), get(obj, 'scan'), get(obj, 'session'), ...
        get(obj, 'run'), pdf_id{ii});
%     [stat, msg] = unix(run_remover);
    unix(run_remover);
end