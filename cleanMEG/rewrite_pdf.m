function rewrite_pdf(newData,labels,source,newPref)
% newData is the data to be written, rows for channels. it can be only
% EEG, only MEG, both, just few channels, doesn't matter. make sure the
% labels match the newData (names of channels in a cell array).
% source is a name of 4D pdf file
% newPref is the new prefix, 'rw' by default
PWD=pwd;
if ~exist('source','var')
    source=[];
end
if isempty(source)
    source=ls('c,*');
    source=source(1:end-1);
    %source='c,rfhp0.1Hz';
end
[pathstr, name, ext] = fileparts(source);
source=[name,ext];


if ~isempty(pathstr)
    cd (pathstr);
end
if ~exist('newPref','var')
    newPref='';
end
if isempty(newPref)
    newFile=['rw_',source];
elseif ~isempty(strfind(source,'_'))
    newFile=[newPref,',',source];
else
    newFile=[newPref,'_',source];
end
if isunix
    eval(['!cp ',source,' ',newFile])
else
    copyfile(source,newFile)
end
pdf=pdf4D(source);
pdf2=pdf4D(newFile);
chi = channel_index(pdf, {'meg' 'ref' 'TRIGGER' 'RESPONSE' 'UACurrent' 'eeg' 'EXTERNAL'}, 'name');
hdr = get(pdf,'header');
lat = [1 hdr.epoch_data{1,1}.pts_in_epoch];
chn = channel_name(pdf, chi);
data = read_data_block(pdf, lat, chi);
if ~exist('labels','var')
    labels=[];
end
if isempty(labels)
    chiIn=channel_index(pdf, 'meg','name');
    labels = channel_label(pdf, chiIn);
end
% for i=1:size(chi,2)
%     dataC(chi(i),:)=data(i,:);
% end
dataC(chi,:)=data; % same but faster
clear data
%% example:
% chi(2)=10;
% chn(2)='A2';
% this means that in data row 2 has channel A2
% in dataC row 10 has channel A2
% chi(276)=307;chn(276)='E2'; on dataC row 307 has E2
%%
newData=single(newData);
for li=1:size(newData,1)
    dataC(chi(find(strcmp(labels(li),chn))),:)=newData(li,:); %#ok<FNDSB>
    %dataC(chi(i),1:size(newData,2))=data(find(strcmp((chn(i)),labels)),:);
end
clear newData
write_data_block(pdf2, dataC, 1);
cd(PWD);
end