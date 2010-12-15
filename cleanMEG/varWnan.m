function m = varWnan(x,dim)
% compute the vatiance even if there are NaNs along dimension dim
% m = varWnan(x,dim);
%
% x   - the data, either array or up to 4D matrix
% dim - dimension along which to compute, if
%       x  is a singleton the var along its greater then 1 dimension is used.
% v   - let N be the number of non NaN
%       elements.  If N is1 or 0 a NaN is returned

% Jan 2008 MA
% UPDATES
% aug-2009 extended for more then 2D  MA

%% initialize
if ~exist('dim','var'), dim=1; end
if isempty(dim), dim=1; end

if isempty(find(isnan(x),1))  % no nan do as usual
    m = var(x,dim);
    return
end

xDim = size(x);
nDim = length(xDim);
if dim>nDim
    error('MATLAB:statistc:mismatchDim','dim cannot be larger then dimensions of x')
end

%% treat one dimensional array
if nDim==2 && xDim(1)==1
    y = x(~isnan(x));
    m = var(y);
    return
elseif nDim==2 && xDim(2)==1
    y = x(~isnan(x));
    m = var(y);
    return
end
switch dim
    case 1
        switch nDim
            case 2
                s=sprintf('m=nan(1,%d);',xDim(2));
            case 3
                s=sprintf('m=nan(1,%d,%d);',xDim(2:end));
            case 4
                s=sprintf('m=nan(1,%d,%d,%d);',xDim(2:end));
        end
    case 2
        switch nDim
            case 2
                s=sprintf('m=nan(%d,1);',xDim(1));
            case 3
                s=sprintf('m=nan(%d,1,%d);',xDim([1,3]));
            case 4
                s=sprintf('m=nan(%d,1,%d,%d);',xDim([1,3:end]));
        end
    case 3
        switch nDim
            case 3
                s=sprintf('m=nan(%d,%d,1);',xDim(1:2));
            case 4
                s=sprintf('m=nan(%d,%d,1,%d);',xDim([1:2,4:end]));
        end
    case 4
        s=sprintf('m=nan(%d,%d,%d);',xDim(1:3));
    otherwise
        error('MATLAB:statistc:inapropriatevalue','dim can be at most 4')
end
eval(s)

%% compute the means
% this is UGLY code
% we know that x has 2 or more dimensions
switch nDim
    case 2
        switch dim
            case 1
                for ii = 1:xDim(2)
                    y = x(:,ii);
                    y(isnan(y))=[];
                    if length(y)>1
                        m(1,ii)=var(y);
                    else
                        m(1,ii)=NaN;
                    end
                end
                return
            case 2
                for ii = 1:xDim(1)
                    y = x(ii,:);
                    y(isnan(y))=[];
                    if length(y)>1
                        m(ii,1)=var(y);
                    else
                        m(ii,1)=NaN;
                    end
                end
                return
        end
    case 3
        switch dim
            case 1
                for ii = 1:xDim(2)
                    for jj = 1:xDim(3);
                        y = x(:,ii,jj);
                        y(isnan(y))=[];
                        if length(y)>1
                            m(1,ii,jj)=var(y);
                        else
                            m(1,ii,jj)=NaN;
                        end
                    end
                end
                return
            case 2
                for ii = 1:xDim(1)
                    for jj = 1:xDim(3);
                        y = x(ii,:,jj);
                        y(isnan(y))=[];
                        if length(y)>1
                            m(ii,1,jj)=var(y);
                        else
                            m(ii,1,jj)=NaN;
                        end
                    end
                end
                return
            case 3
                for ii = 1:xDim(1)
                    for jj = 1:xDim(2);
                        y = x(ii,jj,:);
                        y(isnan(y))=[];
                        if length(y)>1
                            m(ii,jj,1)=var(y);
                        else
                            m(ii,jj,1)=NaN;
                        end
                    end
                end
                return
        end
    case 4
        switch dim
            case 1
                for ii = 1:xDim(2)
                    for jj = 1:xDim(3);
                        for kk = 1:xDim(4)
                            y = x(:,ii,jj,kk);
                            y(isnan(y))=[];
                            if length(y)>1
                                m(1,ii,jj,kk)=var(y);
                            else
                                m(1,ii,jj,kk)=NaN;
                            end
                        end
                    end
                    return
                end
            case 2
                for ii = 1:xDim(1)
                    for jj = 1:xDim(3);
                        for kk = 1:xDim(4)
                            y = x(ii,:,jj,kk);
                            y(isnan(y))=[];
                            if length(y)>1
                                m(ii,1,jj,kk)=var(y);
                            else
                                m(ii,1,jj,kk)=NaN;
                            end
                        end
                    end
                    return
                end
            case 3
                for ii = 1:xDim(1)
                    for jj = 1:xDim(2);
                        for kk = 1:xDim(4)
                            y = x(ii,jj,:,kk);
                            y(isnan(y))=[];
                            if length(y)>1
                                m(ii,jj,1,kk)=var(y);
                            else
                                m(ii,jj,1,kk)=NaN;
                            end
                        end
                    end
                    return
                end
            case 4
                for ii = 1:xDim(1)
                    for jj = 1:xDim(2);
                        for kk = 1:xDim(3)
                            y = x(ii,jj,kk,:);
                            y(isnan(y))=[];
                            if length(y)>1
                                m(ii,jj,kk)=var(y);
                            else
                                m(ii,jj,kk)=NaN;
                            end
                        end
                    end
                    return
                end
        end
end
