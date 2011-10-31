function plotOneFile (DD)
%
%  plot DD for one spike
% DD       - a array of structures each for one single unit

% Aug-2010 generated from plotOneDD
% UPDATES
% 

%% initialize
% colorOrder='rgbmcky';
% find the time range
numRasters = size(DD , 2);
maxT = 0;
for ii = 1:numRasters
    if max(DD(ii).raster.times)>maxT
        maxT = max(DD(ii).raster.times);
    end
end
maxT = 1.05*maxT; % add alittle space at the end
bias = 0;

% % get the table of events
% evList=[];
% jj=1;
% nn=1;
% for ii=1:length(DD(1).raster)
%     evnts = DD(1).raster(ii).evnts;
%     for kk=1:size(evnts ,1)
%         allEvnts(nn,1:2) = evnts(kk,:);
%         allEvnts(nn,3) = ii;  % raster Number
%         nn=nn+1;
%         if ~isnan(evnts(kk,1))
%             evList(jj) = evnts(kk,1);
%             jj = jj+1;
%         end
%     end
% end
% evList = unique(evList);
% 
% % if needed sort the rasters by event times
% if ~isempty(sortBy)
%     I = find(allEvnts(:,1)==sortBy);
%     [sortedTimes, sortingIndx] = sort (allEvnts(I,2));
%     orderWith = allEvnts(I(sortingIndx),3);
%     order = orderWith;
%     jj = length(order)+1;
%     for ii = 1:numRasters
%         if isempty(find(orderWith==ii)), order(jj) = ii; jj=jj+1; end
%     end
% else
%     order = [1:numRasters];
% end
% 
% %
% if spikes
    for ii = 1:numRasters
        times = DD(ii).raster.times;
        if ~isnan(times(1))
            [x,y] = prepareOneRaster(times);
            y = 0.8*y + bias;
            plot(x,y,'LineWidth', 1, 'color', 'k');
            hold on
            text(-1,(y(1)+y(2))/2,DD(ii).name)
        end
        bias = bias+1;
    end
%     ylabel(num2str(DD(1).name));
% else % plot the events
%     for ii = 1:numRasters
%         thisEvnt = DD(1).raster(order(ii)).evnts;
%         if ~isnan(thisEvnt(1))
%             for jj = 1:size(thisEvnt,1)
%                 [x,y] = prepareOneRaster(thisEvnt(jj,2));
%                 y = 1.2*y + bias;
%                 colorIndx = mod(find(evList==thisEvnt(jj,1))-1,length(colorOrder))+1;
%                 colorC = colorOrder(colorIndx);
%                 plot(x,y,'LineWidth', 2, 'color', colorC);
%                 plot(x,y+1,'LineWidth', 2, 'color', colorC);
%                 plot(x,y+2,'LineWidth', 2, 'color', colorC);
%             end
%         end
%         bias = bias+1;
%     end
% end

%% wrap up
set(gca,'TickDir' , 'out')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);

set(gca,'YLim' , [-0.2,numRasters])
set(gca,'XLim',[0,maxT])
title(DD(1).title)
xlabel('Time,  s.')

return
