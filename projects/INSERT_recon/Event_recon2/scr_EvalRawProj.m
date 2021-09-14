close all;
clear all;

path = '/media/ashley/My Passport/TestLRF/PERA_PlanarReconstructionAlgorithm/PeraScripts/projects/INSERT_recon/Event_recon2/Database_Reconstructions/DOI/AllEvents/UniEW/20190312/EW10Per/Axial';

folder = uigetdir(path);
files = dir(fullfile(folder,'*.mat'));
NumFiles = length(files);
data = cell(1,NumFiles);
% Load all data from folder



%%
for i = 1:NumFiles
    
    filename = files(i).name;
    filepath = [files(i).folder,'/'];
    data{i} = open([filepath, filename]);
    try
    data{i}.NodeData = data{i}.NodeData(~cellfun('isempty',data{i}.NodeData));
    data{i}.NodeData = reshape(data{i}.NodeData, [20,4]);
    catch
        disp(['check file' int2str(i)]);
    end
end

% %%
% for i= 1:NumFiles
%     data{i}.NodeData = reshape(data{i}.NodeData, [20,4]);
% end

%%
doi = input('DOI data? Y (1) or N (0)');
NumNodes = length(data{1}.NodeData);
if NumNodes == 4 && doi == 1
    NumNodes = 1;
end
% Single Node Data
if NumNodes == 1
    if doi == 1
        figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:NumFiles
            TotEven = zeros(4,1);
            
            for j = 1:4
                TotEven(j) =sum(sum(data{i}.NodeData{j}));
            end
            subplot(4,5,i)
            plot(1:4, TotEven);
            title(int2str(filenum));
                       
        end
                
        check = 1;
        while check == 1
            
            %filenum = input('Please pick a file:');
            for filenum = 1:NumFiles

                figure('units','normalized','outerposition',[0 0 1 1]);
                for j =1:4
                    subplot(2,2,j);
                    imagesc(data{filenum}.NodeData{j});
                end
            end
            check = input('Check another node? Y(1) N(0)');
        end
    else
        figure('units','normalized','outerposition',[0 0 1 1]);
        TotEven = zeros(NumFiles,1);
        for i = 1:NumFiles
            TotEven(i) =sum(sum(data{i}.NodeData));
        end
        plot(1:NumFiles, TotEven);
        
        check = 1;
        while check == 1
            
            %filenum = input('Please pick a file:');
            for filenum = 1:NumFiles
            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(data{filenum}.NodeData);
            end
            check = input('Check another node? Y(1) N(0)');
        end
        
    end
         if doi == 1
            figure('units','normalized','outerposition',[0 0 1 1]);
            for filenum =1:NumFiles
            imglay = data{filenum}.NodeData;
            img = zeros(258,506);
            for j = 1:4
                img = img + imglay{j}; 
            end
            subplot(4,5,filenum);
            imagesc(img);
            title(int2str(filenum));
            end
            
        else
            figure('units','normalized','outerposition',[0 0 1 1]);
            for i =1:20
            img = data{filenum}.NodeData(i,:);
            subplot(4,5,i);
            imagesc(sum(img{i},3));
            title(int2str(i)); 
            end
            
        end
else
    contin = 1;
    while contin == 1
    filenum = input('Please pick a file number:');
    
    if doi == 1
        figure('units','normalized','outerposition',[0 0 1 1]);
        for j = 1:NumNodes
            NodeData = data{filenum}.NodeData(j,:);
            layer = ones(1,4);
            for k = 1:4
                layer(k) = sum(sum(NodeData{k}));
            end
            subplot(4,5,j);
            plot(1:4,layer);
            title(int2str(j))
            
        end        
        
    else
        NodeData = data{filenum}.NodeData{1};
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        x = 1:NumNodes;
        plot(x,sum(sum(NodeData)));
        
    end
    
    
    
    check = 1;
    while check == 1
        
        %node = input('Please pick a node:');
        for node = 1:20
        if doi == 1
            figure('units','normalized','outerposition',[0 0 1 1]);
            
            imglayer = data{filenum}.NodeData(node,:);
            for j =1:4
                subplot(2,2,j);
                imagesc(imglayer{j});
                xlabel(int2str(j));
            end
           title(int2str(node))
        else
            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(data{filenum}.NodeData{node}(:,:));
        end
        end
        check = input('Check another node? Y(1) N(0)');
    end
    
    contin = input('Continue to another file? Y(1) N(0)');
        if doi == 1
            figure('units','normalized','outerposition',[0 0 1 1]);
            for i =1:20
            imglay = data{filenum}.NodeData(i,:);
            img = zeros(258,506);
            for j = 1:4
                img = img + imglay{j}; 
            end
            subplot(4,5,i);
            imagesc(img);
            title(int2str(i)); 
            end
%            title('total data')
        else
            figure('units','normalized','outerposition',[0 0 1 1]);
            for i =1:20
            img = data{filenum}.NodeData(i,:);
            subplot(4,5,i);
            imagesc(sum(img{i},3));
            title(int2str(i)); 
            end
            
        end
    end
    %%

    
end
