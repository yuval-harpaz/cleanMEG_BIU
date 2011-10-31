function dsel_callback(obj, ev)

% h = guidata(obj);
% data_path = getappdata(h.fig, 'data_path');

switch get(obj, 'Tag')
    case 'init'
        %call from dsel
        update_list(obj, 2);
    case 'add_path'
        h = guidata(obj);
        %get current stage
        data_path = get_data_path;
        %there is no stage use current directory
        if isempty(data_path)
            data_path = pwd;
        end
        %prompt for new data path
        data_path = uigetdir(data_path, 'Select Data Directory');
        %cancel
        if ~ischar(data_path)
            return
        end

        %get current lists
        list = get(h.block(1).block_list, 'String');
        list_bkp = getappdata(h.block(1).block_list, 'List');
        select = getappdata(h.block(1).block_list, 'Select');

        %do not add if already there
        list_match = strcmp(list_bkp, data_path);
        if any(list_match)
            msgbox(['Path ' data_path ' is already in the list'], ...
                'Data Directory');
            return
        end
        
        %add item to the list
        list{end+1} = sprintf('+%s', data_path);

        %update backup list
        list_bkp{end+1} = data_path;

        %store new item in "select" vector
        select(end+1) = true;

        %store new lists
        set(h.block(1).block_list, ...
            'String', list, ...
            'Value', length(list));
        setappdata(h.block(1).block_list, 'List', list_bkp);
        setappdata(h.block(1).block_list, 'Select', select);
        setappdata(h.block(1).block_list, 'Path', list_bkp(select));
        
        %save data_path and save_path in config file
        save_config(obj);
        
        %update patient list
        %(this will also update scan, session, run and pdf lists)
        update_list(obj, 2);
    case 'remove_path'
        %remove current item from list
        h = guidata(obj);
        list = get(h.block(1).block_list, 'String');
        if isempty(list)
            return
        end
        value = get(h.block(1).block_list, 'Value');
        list_bkp = getappdata(h.block(1).block_list, 'List');
        select = getappdata(h.block(1).block_list, 'Select');
        list(value) = [];
        list_bkp(value) = [];
        select(value) = [];
        set(h.block(1).block_list, ...
            'String', list, ...
            'Value', 1);
        setappdata(h.block(1).block_list, 'List', list_bkp);
        setappdata(h.block(1).block_list, 'Select', select);        
        setappdata(h.block(1).block_list, 'Path', list_bkp(select));    

        %save data_path and save_path in config file
        save_config(obj);

        %update patient list
        update_list(obj, 2);
    case 'tool_list'
        h = guidata(obj);
        if get(obj, 'Value') == 1
            set(h.remove_tool, 'Enable', 'off');
            set(h.apply, 'String', 'Save List');
        else
            set(h.remove_tool, 'Enable', 'on');
            set(h.apply, 'String', 'Apply Tool');
        end
    case 'add_tool'
        h = guidata(obj);
        tool_list = get(h.tool_list, 'String');
        new_tool = inputdlg('Enter Function Name', 'Add Tool');
        if isempty(new_tool) || isempty(new_tool{1})
            return
        end
        tool_list{end+1} = new_tool{1};
        set(h.tool_list, 'String', tool_list);
        set(h.tool_list, 'Value', numel(tool_list));
        set(h.remove_tool, 'Enable', 'on');
        set(h.apply, 'String', 'Apply Tool');
        save_config(obj);
    case 'remove_tool'
        h = guidata(obj);
        tool_select = get(h.tool_list, 'Value');
        %if Value is 1, "remove tool" should be disable
        %but just in case:
        if tool_select == 1
            return
        end
        tool_list = get(h.tool_list, 'String');
        tool_list(tool_select) = [];
        set(h.tool_list, 'String', tool_list);
        set(h.tool_list, 'Value', 1);
        set(h.apply, 'String', 'Save List');
        set(obj, 'Enable', 'off');
        save_config(obj);
    case 'apply'
        h = guidata(obj);
        %list of selected files
        list = getappdata(h.block(6).block_list, 'Path');
        tool_select = get(h.tool_list, 'Value');
        if tool_select == 1
            if isempty(list)
                errordlg('No Files to Save', 'Save List');
                return
            end
            save_path = getappdata(h.fig, 'Save_Path');
            old_path = pwd;
            cd(save_path)
            [f, p] = uiputfile( ...
                {'*.*', 'All Files (*.*)'}, ...
                'Save File List as');
            cd(old_path);

            if ~ischar(f)
                return
            end

            %if new save_path store as appdata and update config file
            if ~strcmp(save_path, p)
                setappdata(h.fig, 'Save_Path', p);
                save_config(obj);
            end

            fid = fopen(fullfile(p,f), 'w');
            if fid == -1
                return
            end
            for ii = 1:numel(list)
                fprintf(fid, '%s\n', list{ii});
            end
            fclose(fid);
        else
            tool_list = get(h.tool_list, 'String');
            tool_select = get(h.tool_list, 'Value');
            feval(tool_list{tool_select}, list)
        end
    case 'exit'
        %delete dsel figure
        h = guidata(obj);
        delete(h.fig);
    case 'block_list'
        list = get(obj, 'String');
        if isempty(list)
            return
        end
        value = get(obj, 'Value');
        level = getappdata(obj, 'Level');
        select = getappdata(obj, 'Select');
        switch list{value}(1)
            case '+'
                list{value}(1) = '-';
                select(value) = false;
            case '-'
                list{value}(1) = '+';
                select(value) = true;
        end
        set(obj, 'String', list);
        setappdata(obj, 'Select', select);
        %update current list
        update_list(obj, level);
    case 'filter'
        h = guidata(obj);
        level = getappdata(obj, 'Level');
        filter = get(h.block(level).filter, 'String');
        if level == 1
            if isempty(filter)
                filter = '.*';
                set(h.block(1).filter, 'String', filter);
            end
            list_bkp = getappdata(h.block(1).block_list, 'List');
            select = getappdata(h.block(1).block_list, 'Select');
            list = get(h.block(1).block_list, 'String');
            rgxp = regexp(list_bkp, filter);
            for ii = 1:numel(rgxp)
                if isempty(rgxp{ii})
                    select(ii) = false;
                    list{ii}(1) = '-';
                else
                    select(ii) = true;
                    list{ii}(1) = '+';
                end
            end
            setappdata(h.block(1).block_list, 'Select', select);
            set(h.block(1).block_list, 'String', list);
%             print_list(obj, 1);
        end
        if isempty(filter)
            filter = '*';
            set(h.block(level).filter, 'String', filter);
        end
        %update current list
        update_list(obj, level);
   case 'select_all'
        h = guidata(obj);
        level = getappdata(obj, 'Level');
        block_list = h.block(level).block_list;
        select = getappdata(block_list, 'Select');
        list = get(block_list, 'String');
        for ii = 1:numel(list)
            list{ii}(1) = '+';
            select(ii) = true;
        end
        set(block_list, 'String', list);
        setappdata(block_list, 'Select', select);
        %update current list
        update_list(obj, level);
    case 'select_none'
        h = guidata(obj);
        level = getappdata(obj, 'Level');
        block_list = h.block(level).block_list;
        select = getappdata(block_list, 'Select');
        list = get(block_list, 'String');
        for ii = 1:numel(list)
            list{ii}(1) = '-';
            select(ii) = false;
        end
        set(block_list, 'String', list);
        setappdata(block_list, 'Select', select);
        %update current list
        update_list(obj, level);
    otherwise
        error('Wrong Tag: %s', get(obj, 'Tag'))
end

%%%%%%%%%%%%%%%     Update List     %%%%%%%%%%%%%%%

function update_list(obj, level)

h = guidata(obj);

if level == 1
    %first update level 1 and then switch to level 2
    list = getappdata(h.block(1).block_list, 'List');
    select = getappdata(h.block(1).block_list, 'Select');
    setappdata(h.block(1).block_list, 'Path', list(select));
    level = 2;
end

%path for level above current
list_above = getappdata(h.block(level-1).block_list, 'Path');

% %debuging:
% print_list(obj, level-1);

[list, full_list] = sub_dir(obj, level, list_above);

%add to the list and get current selections
sel = replace_list(obj, level, list);

% %debuging:
% print_list(obj, level);

%store current full list
full_list = full_list(sel);
setappdata(h.block(level).block_list, 'Path', full_list);

if level == 6
    set(h.fig, 'Name', sprintf('%s : %d files', ...
        getappdata(h.fig, 'Name'), numel(full_list)));
end

% %debuging:
% print_list(obj, level);

if level < numel(h.block)
    update_list(obj, level+1)
end

%%%%%%%%%%%%%%%     Replace List     %%%%%%%%%%%%%%%

function sel = replace_list(obj, level, new_list)

%try to add unique items to the list
h = guidata(obj);

%select for the new list
sel = true(1,numel(new_list));

%current list
old_list = getappdata(h.block(level).block_list, 'List');
old_select = getappdata(h.block(level).block_list, 'Select');

%unselected items
old_unlist = old_list(~old_select);

for ii = 1:numel(new_list)
    %check if it was selected in pld list
    list_match = strcmp(old_unlist, new_list{ii});
    if any(list_match)
        %this item was unselected
        sel(ii) = false;
    end
end

%remove repetitive items
[new_list, list_ind] = unique(new_list);

select = sel(list_ind);

list = {};
for ii = 1:numel(new_list)
    if select(ii)
        %add selected item to the list
        list{ii} = sprintf('+%s', new_list{ii});
    else
        %add unselected item to the list
        list{ii} = sprintf('-%s', new_list{ii});
    end

end

%store new lists
%try to keep same listbox value (same row), otherwise set to first or last
old_value = get(h.block(level).block_list, 'Value');
set(h.block(level).block_list, 'String', list);
if old_value == 0 && numel(list)>0
    set(h.block(level).block_list, 'Value', 1);
end
if old_value > numel(list)
    set(h.block(level).block_list, 'Value', numel(list));
end

setappdata(h.block(level).block_list, 'List', new_list);
setappdata(h.block(level).block_list, 'Select', select);

%%%%%%%%%%%%%%%     Sub Dir     %%%%%%%%%%%%%%%

function [list, full_list] = sub_dir(obj, level, list_above)

list = {};
full_list = {};

h = guidata(obj);

%current level filter
filter = get(h.block(level).filter, 'String');

no = 0;
for ii = 1:numel(list_above)
    list_dir = dir(fullfile(list_above{ii}, filter));
%     list_dir = dir([list_above{ii} filesep filter]);
    if isempty(list_dir)
        continue
    end
    for di = 1:numel(list_dir)
        if list_dir(di).isdir && level == 6
            continue
        end
        if list_dir(di).isdir && list_dir(di).name(1) =='.'
            continue
        end
        if ~list_dir(di).isdir && level < 6
            continue
        end
            no = no + 1;
            list{no} = list_dir(di).name;
            full_list{no} = fullfile(list_above{ii}, list_dir(di).name);
    end
end

%%%%%%%%%%%%%%%     Save Dsel Config     %%%%%%%%%%%%%%%

function save_config(obj)

h = guidata(obj);

config_name  = getappdata(h.fig, 'Config');

%this is PC:
if isempty(config_name)
    return
end

data_path = getappdata(h.block(1).block_list, 'List');
save_path = getappdata(h.fig, 'Save_Path');
tool_list = get(h.tool_list, 'String');
config = struct(...
    'data_path', {data_path}, ...
    'tool_list', {tool_list}, ...
    'save_path', save_path);
save(config_name, 'config', '-mat');

%%%%%%%%%%%%%%%     Print List (debuging)     %%%%%%%%%%%%%%%

function print_list(obj, level)

h = guidata(obj);

path_list = getappdata(h.block(level).block_list, 'Path');
list_bkp = getappdata(h.block(level).block_list, 'List');
select = getappdata(h.block(level).block_list, 'Select');
list = get(h.block(level).block_list, 'String');

fprintf('====================================================\n');
fprintf('Level: %d, List size: %d(%d:%d), Path size: %d\n', ...
    level, numel(list), numel(list_bkp), sum(select), numel(path_list));
for ii = 1:numel(list)
    fprintf('list: %d: %s : %s\n', select(ii), list{ii}, list_bkp{ii});
end
for ii = 1:numel(path_list)
    fprintf('path: %s\n', path_list{ii});
end