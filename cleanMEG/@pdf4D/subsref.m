function [B, varargout] = subsref(A,S)

% if length(S)>1
%     error('Wrong reference for class "pdf4D"';
% end

%to handle "Too many output arguments." error in case of A is an array of
%pdf4D objects
varargout = cell(size(A));

B = A;
for n = 1:length(S)

    switch S(n).type
        
        case '.'
            
            switch lower(S(n).subs)
                
                case 'header'
                    %get will get header from structure or from file
                    B = get(B, 'header');
                    
                case 'config'
                    B = get(B, 'config');
                    
                case {'hs' 'headshape'}
                    B = get(B, 'headshape');
                    
                otherwise
                    error('Wrong field name "%s"', S(n).subs)
            end
            
        case '()'
            B = B(S(n).subs{1});
            
        otherwise
            error('Wrong subs type "%s"', S(n).type)
    end
end