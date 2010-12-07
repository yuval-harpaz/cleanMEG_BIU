function pdf = duplicate(source, destination)

% DUPLICATE pdf = duplicate(source, destination)
% source - pdf to be duplicated
% destination - if string then new pdf would be created in the same
%       directory as source, if pdf4D object then source would be
%       copied into directory of destination
% pdf_id - filename for the new pdf (the rest of id(s) taken from the obj)
% pdf - pdf4D object for the new pdf

%create dummy pdf (create database entry)
%(just in case write 3 point and not 1)
if isa(destination, 'pdf4D')
    pdf_id = get(source, 'pdf');
    write_pdf(destination, zeros(1,3), pdf_id, {'A1'}, 1, 0);
    %new filename (it obviously works for destination)
    new_pdf = fullfile(fileparts(destination.FileName), pdf_id);
else
    pdf_id = get(destination, 'pdf');
    write_pdf(source, zeros(1,3), pdf_id, {'A1'}, 1, 0);
    %new filename
    new_pdf = fullfile(fileparts(source.FileName), pdf_id);
end

%duplicate current pdf
copyfile(source.FileName, new_pdf, 'f');

%create new pdf4D object
pdf = pdf4D(new_pdf);