function [ eventData ] = FcnNamecPaolo( N_packets_fixed, C_FILE_INPUT, pacchetti_fissi )
%FCNNAMEPAOLO Summary of this function goes here
%   Detailed explanation goes here

% CONSTANTS
% Runtime options
C_DEBUG             = 0; % C_DEBUG mode is not implemented in this script
C_USE_VERIFICATION  = 0; % When test data is acquired we can do verification for the file because the result is determined
C_SHOW_RESULTS      = 1; % Show results in figures
C_USE_PLOT          = 0; % Plot channels in one figure (from 1 to 36)
C_USE_MAP           = 0; % Plot channels in one image (from 1 to 36 - 6x6)
C_PAUSE_TIME_BETWEEN_FRAMES = 0.1; % Pause time between plots in seconds

% Ethernet packet constants
C_ETHERNET_PACKET_IPV_BYTES        = 42;         % Ipv4 header bytes
C_ETHERNET_PACKET_NAMEC_BYTES      = 17;         % NAMEC bytes
C_ETHERNET_PACKET_FIRMWARE_TYPE    = 4;          % FW bytes
C_ETHERNET_PACKET_DATA_HEADER      = 2;          % Data header bytes
    
% SPECT Insert packet constants
C_NAMEC_SIDW_TYPE_HEADER_WORDS     = 4;          % 32 bit words of header
C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS  = 36;         % 16 bit words of datas
C_SORT_DSR0 = [20  2 19  1 22  4 21  3 24  6 23  5 26  8 25  7 28 10 27  9 30 12 29 11 32 14 31 13 34 16 33 15 36 18 35 17]; % Reorder event data
C_SIDW_PACKET_NUMBER               = 93;         % Used only for C_DEBUG the acquisition errors
 
% Initializations
eventTimestamp = uint64(0);


%% Read from file (simulation file or real acquisition)


%% FILE OPERATIONS
display('---> FILE READING');
Number_Data_Packet = (2*C_NAMEC_SIDW_TYPE_HEADER_WORDS+C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS)*C_SIDW_PACKET_NUMBER + C_ETHERNET_PACKET_DATA_HEADER*2;
fid = fopen(C_FILE_INPUT,'r','b');
 
File = fread(fid,'uint16=>uint16','b');
if pacchetti_fissi
    N_packets = N_packets_fixed;
else
    N_packets = length( File)/Number_Data_Packet;
end
% N_packets = 5000;
Step_Data_Header = (2*C_NAMEC_SIDW_TYPE_HEADER_WORDS+C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS)*C_SIDW_PACKET_NUMBER + C_ETHERNET_PACKET_DATA_HEADER*2;



%% Blocco Data header
Index_Data_Header = [];
for index_0 = 1 : 2*C_ETHERNET_PACKET_DATA_HEADER
    Index_Data_Header = [Index_Data_Header,[index_0:Step_Data_Header:N_packets*Step_Data_Header]];
end
    
File_Data_Header = File(Index_Data_Header);
% File_Data_Header = reshape( File_Data_Header, N_packets, 2*C_ETHERNET_PACKET_DATA_HEADER);
% 
% for index_1 = 1 : C_ETHERNET_PACKET_DATA_HEADER
%     j1 = (index_1-1)*2+1;
%     j2 = (index_1-1)*2+2;
%     Packet_Data_Header(:,index_1) = bin2dec([dec2bin(File_Data_Header(:,j1),16),dec2bin(File_Data_Header(:,j2),16)]);
% end
% 
% A1 = repmat(uint32(hex2dec('000000ff')),N_packets,1);
% A2 = repmat(uint32(hex2dec('0000ff00')),N_packets,1);
% A3 = repmat(uint32(hex2dec('00ff0000')),N_packets,1);
% A4 = repmat(uint32(hex2dec('ff000000')),N_packets,1);
% % Event type
% packetDataType(:,4) = char(uint8(bitshift(bitand(double(A1),Packet_Data_Header(:,2)),-0)));
% packetDataType(:,3) = char(uint8(bitshift(bitand(double(A2),Packet_Data_Header(:,2)),-8)));
% packetDataType(:,2) = char(uint8(bitshift(bitand(double(A3),Packet_Data_Header(:,2)),-16)));
% packetDataType(:,1) = char(uint8(bitshift(bitand(double(A4),Packet_Data_Header(:,2)),-24)));

%% Blocco SIDW header
Index_SIDW_Header = [];
for index_packet = 1 : N_packets
    if rem(index_packet,1000)==0
        index_packet
    end
    Counter_packet = (index_packet-1)*Number_Data_Packet;  
    Index_Start = 2*C_ETHERNET_PACKET_DATA_HEADER+1+Counter_packet;
    Step_SIDW_Header = C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS+2*C_NAMEC_SIDW_TYPE_HEADER_WORDS;
    Index_End = Index_Start+Step_SIDW_Header*C_SIDW_PACKET_NUMBER-1;
    Index_SIDW_Header = [Index_SIDW_Header,[Index_Start : Step_SIDW_Header : Index_End]];
end      
 
Index_SIDW_Header_0 = Index_SIDW_Header;
Index_SIDW_Header = [];
for index_0 = 0 : 2*C_NAMEC_SIDW_TYPE_HEADER_WORDS-1
    Index_SIDW_Header = [Index_SIDW_Header, Index_SIDW_Header_0+index_0];
end
clear Index_SIDW_Header_0

Packet_Data_SIDW_Header = File(Index_SIDW_Header);
Packet_Data_SIDW_Header = reshape( Packet_Data_SIDW_Header, C_SIDW_PACKET_NUMBER*N_packets, 2*C_NAMEC_SIDW_TYPE_HEADER_WORDS);
% 
% for index_1 = 1 : C_NAMEC_SIDW_TYPE_HEADER_WORDS
%     index_1
%     j1 = (index_1-1)*2+1;
%     j2 = (index_1-1)*2+2;
%     A = [dec2bin(Packet_Data_SIDW_Header(:,j1),16),dec2bin(Packet_Data_SIDW_Header(:,j2),16)];
%     Packet_Data_SIDW_Header_32bit(:,index_1) = bin2dec(A);
% end
% 
% A1 = repmat(uint32(hex2dec('000000ff')),C_SIDW_PACKET_NUMBER*N_packets,1);
% A2 = repmat(uint32(hex2dec('0000ff00')),C_SIDW_PACKET_NUMBER*N_packets,1);
% A3 = repmat(uint32(hex2dec('00ff0000')),C_SIDW_PACKET_NUMBER*N_packets,1);
% A4 = repmat(uint32(hex2dec('ff000000')),C_SIDW_PACKET_NUMBER*N_packets,1);
% % Event type
% eventDataTypeTemp(:,4) = char(uint8(bitshift(bitand(double(A1),Packet_Data_SIDW_Header_32bit(:,2)),-0)));
% eventDataTypeTemp(:,3) = char(uint8(bitshift(bitand(double(A2),Packet_Data_SIDW_Header_32bit(:,2)),-8)));
% eventDataTypeTemp(:,2) = char(uint8(bitshift(bitand(double(A3),Packet_Data_SIDW_Header_32bit(:,2)),-16)));
% eventDataTypeTemp(:,1) = char(uint8(bitshift(bitand(double(A4),Packet_Data_SIDW_Header_32bit(:,2)),-24)));
% % Event timestamp
% eventTimestampMsbTemp = bitand(double(uint32(hex2dec('0000ffff'))),Packet_Data_SIDW_Header_32bit(:,3));



%% Blocco dati
Index_SIDW_Data = [];

for index_packet = 1 : N_packets
    if rem(index_packet,1000)==0
        index_packet
    end
    Counter_packet = (index_packet-1)*Number_Data_Packet;  
    Index_Start = 2*C_ETHERNET_PACKET_DATA_HEADER+2*C_NAMEC_SIDW_TYPE_HEADER_WORDS+1+Counter_packet;
    Index_End = Index_Start+Step_SIDW_Header*C_SIDW_PACKET_NUMBER-1;
    Index_SIDW_Data = [Index_SIDW_Data,[Index_Start : Step_SIDW_Header : Index_End]];
end      
 
Index_SIDW_Data_0 = Index_SIDW_Data;
Index_SIDW_Data = [];
for index_0 = 0 : C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS-1
    Index_SIDW_Data = [Index_SIDW_Data, Index_SIDW_Data_0+index_0];
end        
clear Index_SIDW_Data_0

Packet_Data = File(Index_SIDW_Data);
Packet_Data = reshape( Packet_Data, C_SIDW_PACKET_NUMBER*N_packets, C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS);

column = C_SORT_DSR0;
eventData(:,column) = Packet_Data(:,1:C_NAMEC_SIDW_DSR0_TYPE_DATA_WORDS);
fclose(fid);


end

