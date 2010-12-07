function A = subsasgn(A,S,B)

if length(S)>1
    error('Wrong assignment for class "pdf4D"');
end

switch S.type
    case '.'
        switch S.subs
            case 'header'
                %'set' might check if B is a header (for now it does not)
                A = set(A, 'header', B);
            case 'config'
                A = set(A, 'config', B);
            case {'hs' 'headshape'}
                A = set(A, 'headshape', B);
            otherwise
                error('Wrong field name "%s"', S.subs)
        end
    otherwise
        error('Wrong subs type "%s"', S.type)
end