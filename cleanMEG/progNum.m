function progNum(ii)
% shows progress (channel number)
if ii==1
    fprintf(num2str(ii))
else
    txt=[repmat('\b',1,length(num2str(ii-1))),num2str(ii)];
    fprintf(txt)
end