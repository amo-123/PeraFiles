function [ii] = MS_getfilesize(File_Name,File_Path,num_events)

disp(strcat('Loading .data File:',{' '},File_Name));

%initializations and parameters
time_res = 25e-9; 
disp(strcat('time resolution =',32,num2str(time_res)));

%% File opening

fid = fopen([File_Path,File_Name],'r','b'); %windows

%check file dimension
fseek(fid,0,'eof');
FileSize=ftell(fid); %location of the file position indicator is indicated in bytes from the beginning of the file
fseek(fid,0,'bof');

file = fread(fid,8,'uint16=>uint16','b');
fseek(fid,0,'bof');

%evaluate if Clinical or Preclinical
if file(2) ==  51 && file(8) == [21041]
    %Clinical
%     modality = 'Clinical';
     pSize = 51;
    channels = 72;
%     text_header = 21041;
elseif file(2) ==  93 && file(8) == [21040]
    %Preclinical
%     modality = 'Preclinical';
     pSize = 93;
    channels = 36;
%     text_header = 21040;
else
%     modality = 'Error';
end
% disp(strcat('modality=',num2str(modality)));

%estimated number of events, considering 16bits words
% estimated_n_events=floor((pSize*FileSize/2)/((8+channels)*pSize+4));

file_pos = 0:2*(num_events*(8+channels)+4):FileSize;

if length(file_pos) < 3
    estimated_n_events = floor((pSize*FileSize/2)/((8+channels)*pSize+4));
    file_pos = [0,2*(estimated_n_events*(8+channels)+4),FileSize];
end

ii = length(file_pos);
