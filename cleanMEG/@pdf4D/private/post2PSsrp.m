function pcell = post2PSsrp(post)

%convert long string into cell array (use '@' as separator)

if isempty(post)
    pcell = {''};
else
    lend = find(post=='@');
    lstart = [1 lend+1];
    lend = [lend - 1 numel(post)];
    for li = 1:numel(lend)
        pcell{li} = post(lstart(li):lend(li));
    end
end
