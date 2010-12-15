function pdf = pdf4D(varargin)

%PDF4D 4D pdf class constructor.
%   p = PDF4D(filename) creates a pdf object for the filename
%
%Distributed under GNU GENERAL PUBLIC LICENSE (see license.txt)
%
%Copyleft 2006, eugene.kronberg@uchsc.edu

switch nargin
    case 0
       % if no input arguments, create a default object
       init_data;
       pdf = class(pdf,'pdf4D');
    case 1
        switch class(varargin{1})
            case 'pdf4D'
                % if single argument of class 'pdf4D', return it
                pdf = varargin{1};
            case 'cell'
                % array of pdf names, return array of pdf4D objects
                names = varargin{1};
                for ind = 1:length(names)
                    pdf(ind) = pdf4D(names{ind});
                end
            case 'char'
                % if single argument of class 'char', this is pdf filename
                if ispdf(varargin{1})
                    init_data;
                    pdf.FileName = varargin{1};
                    id = path2PSsrp(varargin{1});
                    if length(id) == 5
                        pdf.ID.Patient = id{1};
                        pdf.ID.Scan = id{2};
                        pdf.ID.Session = id{3};
                        pdf.ID.Run = id{4};
                        pdf.ID.PDF = id{5};
                    end
                    pdfdir = fileparts(varargin{1});
                    config = fullfile(pdfdir, 'config');
                    if isconfig(config)
                        pdf.ConfigName = config;
                    end
                    hs_file = fullfile(pdfdir, 'hs_file');
                    if ishs(hs_file)
                        pdf.HeadShapeName = hs_file;
                    end
                    pdf = class(pdf,'pdf4D');
                else
                      %uncomment for testing only:
%                     init_data;
%                     pdf.FileName = varargin{1};
%                     pdf = class(pdf,'pdf4D');
                    error('Bad pdf name')
                end
            otherwise
                  error('Wrong argument type')
        end 
    otherwise
       error('Wrong number of input arguments')
end

    function init_data
       obj_path = fileparts(mfilename('fullpath'));
       fid=fopen(fullfile(obj_path, 'version.txt'));
       pdf.Version = fgetl(fid);
       fclose(fid);
       pdf.ID = struct( ...
           'Patient', '', ...
           'Scan', '', ...
           'Session', '', ...
           'Run', '', ...
           'PDF', '');
       pdf.FileName = '';
       pdf.Header = [];
       pdf.ConfigName = '';
       pdf.Config = [];
       pdf.HeadShapeName = '';
       pdf.HeadShape = [];
    end
end