function val = get(obj,propName)

% GET    Get pdf4D object properties.
%    V = GET(OBJ,'PropertyName') returns the value of the specified
%    property for the pdf4D object OBJ.
%  
%    GET(OBJ) displays all property names for the pdf4D object.

if nargin == 1
    fprintf([ ...
        '\tver | version : pdf4D object version\n' ...
        '\tpdffound | filefound : true | false\n' ...
        '\tconfigfound : true | false\n' ...
        '\thsfound | headshapefound : true | false\n' ...
        '\tpdfname | filename : full pdf filename\n' ...
        '\tconfigname : full config filename\n' ...
        '\thsname | headshapename : full headshape filename\n' ...
        '\tid | ids : id structure\n' ...
        '\tpatient | patientid : patient id\n' ...
        '\tscan | scanid : scan id\n' ...
        '\tsession | sessionid : session id\n' ...
        '\trun | runid : run id\n' ...
        '\tpdf | pdfid : pdf id\n' ...
        '\theader : pdf header (structure)\n' ...
        '\tconfig : config info (structure)\n' ...
        '\ths | headeshape : head shape info (structure)\n' ...
        '\thsindex : index points from hs_file (structure)\n' ...
        '\thspoint | hspoints : head shape points from hs_file\n' ...
        '\tdf | dataformat : data format (string)\n' ...
        '\tsp | sampleperiod : sample period [sec]\n' ...
        '\tdr | digitizationrate : digitization rate [Hz]\n' ...
        '\ttrigger : trigger event\n' ...
        '\tfirst_latency : trigger start latency\n' ...
        '\tepoch_duration : epoch duration\n' ...
        '\tpts_in_epoch : total number of points in one epoch\n' ...
        '\ttotal_epochs : total number of epochs in pdf\n' ...
        '\ttotal_chans | total_channels : total number of channels in pdf\n' ...
        '\n'
        ]);
    return
end

%if input is array of pdf4D objects return cell array of props
if length(obj) > 1
    val = cell(1,length(obj));
    for ind = 1:length(obj)
        val{ind} = get(obj(ind), propName);
    end
    return
end

switch lower(propName)
    case {'ver' 'version'}
        val = obj.Version;
    case {'pdffound' 'filefound'}
        val = ~isempty(obj.FileName);
    case 'configfound'
        val = ~isempty(obj.ConfigName);
    case {'hsfound' 'headshapefound'}
        val = ~isempty(obj.HeadShapeName);
    case {'pdfname' 'filename'}
        val = obj.FileName;
    case 'configname'
        val = obj.ConfigName;
    case {'hsname' 'headshapename'}
        val = obj.HeadShapeName;
    case {'id' 'ids'}
        val = obj.ID;
    case {'patientid' 'patient'}
        val = obj.ID.Patient;
    case {'scan_id' 'scan'}
        val = obj.ID.Scan;
    case {'sessionid' 'session'}
        val = obj.ID.Session;
    case 'session_dir'
        val = session2dir(get(obj, 'sessionid'));
    case {'runid' 'run'}
        val = obj.ID.Run;
    case {'pdfid' 'pdf'}
        val = obj.ID.PDF;
    case 'header'
        if isempty(obj.Header)
           if get(obj, 'FileFound')
                val = read_header(obj.FileName);
           else
                val = [];
           end
        else
            val = obj.Header;
        end
    case 'config'
        if isempty(obj.Config)
           if get(obj, 'ConfigFound')
                val = read_config(obj.ConfigName);
           else
                val = [];
           end
        else
            val = obj.Config;
        end
    case {'hs' 'headshape'}
        if isempty(obj.HeadShape)
           if get(obj, 'HeadShapeFound')
                val = read_hs_file(obj.HeadShapeName);
           else
                val = [];
           end
        else
            val = obj.HeadShape;
        end
    case 'hsindex'
        if isempty(obj.HeadShape)
           if get(obj, 'HeadShapeFound')
                hs = read_hs_file(obj.HeadShapeName);
                val = hs.index;
           else
                val = [];
           end
        else
            val = obj.HeadShape.index;
        end
    case {'hspoint' 'hspoints'}
        if isempty(obj.HeadShape)
           if get(obj, 'HeadShapeFound')
                hs = read_hs_file(obj.HeadShapeName);
                val = hs.point;
           else
                val = [];
           end
        else
            val = obj.HeadShape.point;
        end
    case {'df' 'dataformat'}
        %BTi Data Formats:
        df = {'short' 'long' 'float' 'double'};
        hdr = get(obj, 'header');
        if isempty(hdr)
            val = [];
        else
            val = df{hdr.header_data.data_format};
        end
    case {'sp' 'sampleperiod'}
        hdr = get(obj, 'header');
        if isempty(hdr)
            val = [];
        else
            val = hdr.header_data.sample_period;
        end
    case {'dr' 'digitizationrate'}
        sp = get(obj, 'sampleperiod');
        if isempty(sp) || sp==0
            val = [];
        else
            val = 1/sp;
        end
    case 'trigger'
        hdr = get(obj, 'header');
        if isempty(hdr)
            val = [];
        else
            val = [];
            for ei = 1:hdr.header_data.total_events
                %search for trigger event
                if strcmp(hdr.event_data{ei}.event_name, 'Trigger')
                    val = hdr.event_data{ei};
                    break
                end
            end
        end
    case 'first_latency'
        trig = get(obj, 'trigger');
        if isempty(trig)
            val = [];
        else
            val = -trig.start_lat;
        end
    case 'epoch_duration'
        hdr = get(obj, 'header');
        val = hdr.epoch_data{1}.epoch_duration;
    case 'pts_in_epoch'
        hdr = get(obj, 'header');
        val = hdr.epoch_data{1}.pts_in_epoch;
    case 'total_epochs'
        hdr = get(obj, 'header');
        val = hdr.header_data.total_epochs;
    case {'total_chans' 'total_channels'}
        hdr = get(obj, 'header');
        val = hdr.header_data.total_chans;
    otherwise
        error( ...
            'There is no ''%s'' property in the ''pdf4D'' class', ...
            propName)
end