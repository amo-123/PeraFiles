% % FilterSpec = '.data';
% % %[filename,filepath] = uigetfile([pwd,'\',FilterSpec], 'Select .data file', 'MultiSelect', 'off');
% for i =3:18
%     
% folder = uigetdir;
% files = dir(fullfile(folder,'*.data'));
% leng = length(files);
% L_msk = 0;
% 
% for j = 1:leng
%     
%     filename = files(j).name;
% filepath = [files(j).folder,'/'];
% Uniform = UniSum(filename,filepath,L_msk);
% output_savepath = strcat(pwd,'/Database/UniformStudy/Raw/Uniform_',int2str(i),'_',int2str(j),'.mat');
% save(output_savepath, 'Uniform');
% end
% 
% 
% end

%% 

folder = uigetdir;
files = dir(fullfile(folder,'*.mat'));
leng = length(files);
AllUni = cell(leng,1);
for j = 1:leng
    
    filename = files(j).name;
filepath = [files(j).folder,'/'];

AllUni{j} = load(strcat(filepath,filename),'Uniform');
end

%% 

% Integral Uniformity 
IntUni = zeros(41,20);
for k = 1:41
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform > 2e+8) = 0;
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform < 5e+6) = 0;
    
    for m= 1:20
        maxCount = max(nonzeros(AllUni{k,1}.Uniform(m,:)));
        minCount = min(nonzeros(AllUni{k,1}.Uniform(m,:)));
        try
            IntUni(k,m) = (maxCount - minCount)/(maxCount + minCount)*100;
        catch 
            IntUni(k,m) = nan;
        end
%         if IntUni(k,m) < 5
%             IntUni(k,m) = nan;
%         end
    end
end

%%

% Coefficient of Variance

COV = zeros(41,20);
for k = 1:41
    
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform > 2e+8) = 0;
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform < 5e+6) = 0;
    for m= 1:20
        meanCount = mean(nonzeros(AllUni{k,1}.Uniform(m,:)));
        stdCount = std(nonzeros(AllUni{k,1}.Uniform(m,:)));
        try
            COV(k,m) = (stdCount/meanCount)*100;
        catch
            COV(k,m) = nan;
        end
        %COV(isnan(COV))=100;
    end
end

%% 

% Integral Uniformity 
IntUni = zeros(41,20);
for k = 1:41
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform > 2e+8) = 0;
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform < 5e+6) = 0;
    
    for m= 1:20
        %maxCount = max(nonzeros(UFOVKill{k,1}.Uniform(m,:)));
        %minCount = min(nonzeros(UFOVKill{k,1}.Uniform(m,:)));
        
        maxCount = max(nonzeros(killUni{k,1}.Uniform(m,:)));
        minCount = min(nonzeros(killUni{k,1}.Uniform(m,:)));
        %maxCount = max(nonzeros(UFOV{k,1}.Uniform(m,:)));
        %minCount = min(nonzeros(UFOV{k,1}.Uniform(m,:)));
        try
            IntUni(k,m) = (maxCount - minCount)/(maxCount + minCount)*100;
        catch 
            IntUni(k,m) = nan;
        end
        if IntUni(k,m) < 1
            IntUni(k,m) = nan;
        end
    end
end

%%

% Coefficient of Variance

COV = zeros(41,20);
for k = 1:41
    
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform > 2e+8) = 0;
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform < 5e+6) = 0;
    for m= 1:20
   %     meanCount = mean(nonzeros(UFOVKill{k,1}.Uniform(m,:)));
    %    stdCount = std(nonzeros(UFOVKill{k,1}.Uniform(m,:)));

        meanCount = mean(nonzeros(killUni{k,1}.Uniform(m,:)));
        stdCount = std(nonzeros(killUni{k,1}.Uniform(m,:)));
   
   %     meanCount = mean(nonzeros(UFOV{k,1}.Uniform(m,:)));
    %    stdCount = std(nonzeros(UFOV{k,1}.Uniform(m,:)));
        try
            COV(k,m) = (stdCount/meanCount)*100;
        catch
            COV(k,m) = nan;
        end
        %COV(isnan(COV))=100;
    end
end


%% 

% Integral Uniformity 
IntUni = zeros(41,1);
for k = 1:41
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform > 2e+8) = 0;
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform < 5e+6) = 0;
    
    for m= 1:20
        Uall(m,:) = AllUni{k,1}.Uniform(m,:);
    end
        maxCount = max(nonzeros(Uall));
        minCount = min(nonzeros(Uall));
        
        %maxCount = max(nonzeros(killUni{k,1}.Uniform(m,:)));
        %minCount = min(nonzeros(killUni{k,1}.Uniform(m,:)));
        %maxCount = max(nonzeros(UFOV{k,1}.Uniform(m,:)));
        %minCount = min(nonzeros(UFOV{k,1}.Uniform(m,:)));
        try
            IntUni(k) = (maxCount - minCount)/(maxCount + minCount)*100;
        catch 
            IntUni(k) = nan;
        end
        if IntUni(k) < 1
            IntUni(k) = nan;
        end
end

%%

% Coefficient of Variance

COV = zeros(41,1);
for k = 1:41
    
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform > 2e+8) = 0;
%     AllUni{k,1}.Uniform(AllUni{k,1}.Uniform < 5e+6) = 0;
    for m= 1:20
            Uall(m,:) = UFOVKill{k,1}.Uniform(m,:);
       
    end
            meanCount = mean(nonzeros(Uall));
        stdCount = std(nonzeros(Uall));

  %      meanCount = mean(nonzeros(killUni{k,1}.Uniform(m,:)));
   %     stdCount = std(nonzeros(killUni{k,1}.Uniform(m,:)));
   
   %     meanCount = mean(nonzeros(UFOV{k,1}.Uniform(m,:)));
    %    stdCount = std(nonzeros(UFOV{k,1}.Uniform(m,:)));
        try
            COV(k) = (stdCount/meanCount)*100;
        catch
            COV(k) = nan;
        end
end