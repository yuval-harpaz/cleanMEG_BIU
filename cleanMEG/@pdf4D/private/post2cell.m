function pcell = post2cell(post)

%convert long string into cell array (use '\n' as separator)

if isempty(post)
    pcell = {};
else
    %new line character
    nl = sprintf('\n');
    lend = find(post==nl);%last character for each line
    lstart = [1 lend+1];%first character for each line
    lstart(end) = [];%there is no line past last new-line character
    lend = lend - 1;%remove new-line characters from result
    for li = 1:numel(lend)
        pcell{li} = post(lstart(li):lend(li));
    end
end
