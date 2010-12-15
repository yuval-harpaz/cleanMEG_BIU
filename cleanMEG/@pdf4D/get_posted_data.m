function [pdf, id] = get_posted_data(obj)

% GET_POSTED_DATA [pdf, id] = get_posted_data(obj)
% returns full file name and id structure for the posted data

pdf = {};
id = {};

%get posted selection
post = get_posted_sel(obj);
if isempty(post)
    return
end

%convert long string into cell array (one cell per posted data set)
post = post2cell(post);

%return one cell for each posted file
for p_i=1:numel(post)
    pdf{p_i} = post2path(post{p_i});
    id{p_i} = post2PSsrp(post{p_i});
end