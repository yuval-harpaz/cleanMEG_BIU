function chanindx = channel_index(obj, name, mode)

% CHANNEL_INDEX chanindx = channel_index(obj, name, mode)
% name could be channel name or label, or 'meg', 'eeg', 'ext', 'ref',
%    or cell array of names (something like {'TRIGGER' 'A1'})
% mode could be 'name' (name is a name or to sort by names),
%    'label' (name is a label or to sort by labels)

%check number of inputs
if nargin < 2
    error('Wrong number of inputs');
end

%default - no sorting (or this is channel label)
if nargin < 3
    mode = 'none';
end

%BTi channel types
MEG = 1;
EEG = 2;
REFERENCE = 3;
EXTERNAL = 4;
TRIGGER = 5;
UTILITY = 6;
DERIVED = 7;
SHORTED = 8;

%in case there are no channels
chanindx = [];

switch class(name)
    case 'cell'
%         chanindx = [];
        for n = 1:length(name)
            chanindx = [chanindx channel_index(obj, name{n}, mode)];
        end
    case 'char'
        header = get(obj, 'header');
        config = get(obj, 'config');
        switch name
            
            case {'meg' 'MEG'}
                meg_ch = 0;
%                 chanindx = 0;
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= MEG
                        continue
                    end
                    meg_ch = meg_ch + 1;
                    ch_index(meg_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{meg_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %drop 'A' and sort names numerically
                    ch_name_index(meg_ch) = ...
                        str2num(config.channel_data{chan_no}.name(2:end));
                end
                if meg_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name_index);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end
                
            case {'eeg' 'EEG'}
                eeg_ch = 0;
%                 chanindx = 0;
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= EEG
                        continue
                    end
                    eeg_ch = eeg_ch + 1;
                    ch_index(eeg_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{eeg_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %drop 'E' and sort names numerically
                    ch_name_index(eeg_ch) = ...
                        str2num(config.channel_data{chan_no}.name(2:end));
                end
                if eeg_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name_index);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end
                
            case {'ref' 'reference' 'REFERENCE'}
                ref_ch = 0;
%                 chanindx = 0;%in case there is no ref channels
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= REFERENCE
                        continue
                    end
                    ref_ch = ref_ch + 1;
                    ch_index(ref_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{ref_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %sort names alphanumerically
                    ch_name{ref_ch} = ...
                        config.channel_data{chan_no}.name;
                end
                if ref_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end
                
            case {'ext' 'external' 'EXTERNAL'}
                ext_ch = 0;
%                 chanindx = 0;%in case there is no ref channels
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= EXTERNAL
                        continue
                    end
                    ext_ch = ext_ch + 1;
                    ch_index(ext_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{ext_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %drop 'X' and sort names numerically
                    ch_name_index(ext_ch) = ...
                        str2num(config.channel_data{chan_no}.name(2:end));
                end
                if ext_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name_index);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end

            case {'trig' 'trigger'} %upper case for 'TRIGGER' channel
                trig_ch = 0;
%                 chanindx = 0;%in case there is no ref channels
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= TRIGGER
                        continue
                    end
                    trig_ch = trig_ch + 1;
                    ch_index(trig_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{trig_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %sort names alphanumerically
                    ch_name{trig_ch} = ...
                        config.channel_data{chan_no}.name;
                end
                if trig_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end

            case {'util' 'utility' 'UTILITY'}
                util_ch = 0;
%                 chanindx = 0;%in case there is no ref channels
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= UTILITY
                        continue
                    end
                    util_ch = util_ch + 1;
                    ch_index(util_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{util_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %sort names alphanumerically
                    ch_name{util_ch} = ...
                        config.channel_data{chan_no}.name;
                end
                if util_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end

            case {'drv' 'derived' 'DERIVED'}
                drv_ch = 0;
%                 chanindx = 0;%in case there is no ref channels
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= DERIVED
                        continue
                    end
                    drv_ch = drv_ch + 1;
                    ch_index(drv_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{drv_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %sort names alphanumerically
                    ch_name{drv_ch} = ...
                        config.channel_data{chan_no}.name;
                end
                if drv_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end

            case {'shrt' 'shorted' 'SHORTED'}
                shrt_ch = 0;
%                 chanindx = 0;%in case there is no ref channels
                for ch = 1:header.header_data.total_chans
                    chan_no = header.channel_data{ch}.chan_no;
                    if config.channel_data{chan_no}.type ~= SHORTED
                        continue
                    end
                    shrt_ch = shrt_ch + 1;
                    ch_index(shrt_ch) = ch;
                    %sort lables alphanumerically
                    ch_label{shrt_ch} = ...
                        header.channel_data{ch}.chan_label;
                    %sort names alphanumerically
                    ch_name{shrt_ch} = ...
                        config.channel_data{chan_no}.name;
                end
                if shrt_ch == 0
                    return
                end
                switch mode
                    case 'name'
                        [dum, indx] = sort(ch_name);
                        chanindx = ch_index(indx);
                    case 'label'
                        [dum, indx] = sort(ch_label);
                        chanindx = ch_index(indx);
                    case 'none'
                        chanindx = ch_index;
                    otherwise
                        error('Wrong sort mode: %s', mode)
                end

            otherwise %this is channel name or label
%                 chanindx = 0; %in case we would not find channel
                for ch = 1:header.header_data.total_chans
                    switch mode
                        case 'name'
                            chan_no = header.channel_data{ch}.chan_no;
                            if strcmp(name, ...
                                    config.channel_data{chan_no}.name)
                                chanindx = ch;
                                return
                            end
                        case {'label' 'none'}
                            if strcmp(name, ...
                                    header.channel_data{ch}.chan_label)
                                chanindx = ch;
                                return
                            end
                        otherwise
                            error('Wrong sort mode: %s', mode)
                    end
                end
        end
    otherwise
        error('Wrong class for name: %s', class(name))
end