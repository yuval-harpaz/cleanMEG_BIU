function test = msi_installed

%unix with msi
if isunix
%     [s, w] = unix('echo $STAGE');
    [s, w] = unix('which msi');
    test = ~s;
else
    test = false;
end