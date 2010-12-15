function test(obj)

%test if pdf still exists
%add to all methods which use pdf

if ~ispdf(obj.FileName)
    error('Bad pdf object')
end