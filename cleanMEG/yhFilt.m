function outData = yhFilt(inData,filtObj)
% filter back and forth using the filter object
%   outData = myFilt(inData,filtObj);
% inData  - matrix of nChannels X mDataPoints
% filtObj - filter object as created by fdesign. ... and design
% outData - the filtered data in the raws

% Oct-2008  MA
%  UPDATES
%   Jun-2009  Number of columns must be bigger then number of rows

%% initialize
[rows,columns] = size(inData);
if rows>columns
    warning('more rows than columns, transform inData?'); % yuval
end
numChannels = size(inData,1);
% numPoints = size(inData,2);
outData = zeros(size(inData));

%% filter
for ii=1:numChannels
    y = filter(filtObj, inData(ii,:));
    outData(ii,:) = fliplr(filter(filtObj,fliplr(y)));
    if numChannels>1000
        if ii==round(ii/100)*100
            display(num2str(ii))
        end
    end
end

return

