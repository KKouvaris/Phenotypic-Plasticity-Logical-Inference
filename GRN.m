%clear command window and workspace
clc; clear;

%0-static,
%1-heterogeneous,
%2-plasticity,
%3-plasticity binary,
%4-plasticity test set,
%5-static with env cue
%6-plastic anti-correlation
%7-plastic anti-correlation with a redundant trait - 0 env cue for that trait
%8-plastic AB- BC- -> AC+ , as 7, change every generation
%9-plastic AB+ BC- -> AC- , as 7, change every generation
%10-plastic AB+ BC+ -> AC+ , as 7, change every generation
%11-plastic AB+ AB- -> AB0 , as 7, change every generation
%12-plastic {AB+, {BC+, BC-} -> BC0} -> AC0 , as 7, change every generation
%13-plastic {A><B><C}->{AB0,AC0,BC0}
%14-plastic {A,B,C}->{AB0,AC0,BC0}
%15-static as 9 {ABC: AB+, BC-, AC-}
%16-plastic as 9, non-linear FF, binary ----------------
%17-plastic as 9, non-linear FF, bipolar binary
%18-plastic, {ABCD: AB+ BC-}, S=20, [min_cue,max_cue]=[-1,1] change every generation

N = 1000;
init_t = 0; %50000;
scenario = 18;

delete(gcp('nocreate'));
myCluster = parcluster('local');
myCluster.NumWorkers = 10;
saveProfile(myCluster);
parpool(myCluster);

if (init_t == 0)
    
    numElements = 200;
    
else
    
    %get directory and determine the number of existing replicates
    folderDirectory = findFolder( scenario, N, init_t );
    
    %determine the number of replicates in the respective scenario
    dirList = dir([folderDirectory '\*.mat']);
    numList = NaN(length(dirList),1);
    
    for id = 1:length(dirList)
        %Get the file name (minus the extension)
        [~, f] = fileparts(dirList(id).name);
        if ~isempty(strfind(f,'Pop'))
            findToken = strfind(f,'_');
            numList(id) = str2double(f(findToken(1)+1:findToken(2)-1));
        end
    end
    numList(isnan(numList)) = [];
    numList = sort(numList);
    numElements = length(numList);
    
end

for k = 1:ceil(numElements/10)
    saveObj = [];
    parfor i=(k-1)*10 + 1: min(k*10,numElements)
        if (init_t ~= 0)
            init_Pop = load([folderDirectory '\replicate_' num2str(numList(i)) '_Pop.mat']);
            saveObj = [saveObj GRN_Fun(scenario, init_t,init_Pop.Pop)];
        else
            saveObj = [saveObj GRN_Fun(scenario, init_t)];
        end
    end
    %save
    for i=1:length(saveObj)
        saveMyVariables(scenario, saveObj(i));
    end
end



