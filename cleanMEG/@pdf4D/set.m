function obj = set(obj,varargin)

% SET Set pdf4D properties and return the updated object

if isempty(varargin)
    fprintf([ ...
        'header : attach header structure to the object\n' ...
        'config : attach config structure to the object\n' ...
        'hs | headshape : attach head shape structure to the object\n' ...
        ]);
    return
end

propertyArgIn = varargin;

while length(propertyArgIn) >= 2
    
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    
    switch lower(prop)
    case 'header'
        obj.Header = val;
    case 'config'
        obj.Config = val;
    case {'hs' 'headshape'}
        obj.HeadShape = val;
    otherwise
        error('Wrong pdf4D properties: %s', prop)
    end
end