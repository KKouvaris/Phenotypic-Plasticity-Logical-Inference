function [ initFit_Collection_Plastic, fitDiff_Collection_Plastic ] = evolvabilityDS_v2( scenario, N, t_max )

folderDirectory = findFolder(scenario, N, t_max);

%determine the number of replicates in the respective scenario
dirList = dir([folderDirectory '\*.mat']);
numList = NaN(length(dirList),1);

for id = 1:length(dirList)
    % Get the file name (minus the extension)
    [~, f] = fileparts(dirList(id).name);
    if ~isempty(strfind(f,'Pop'))
        findToken = strfind(f,'_');
        numList(id) = str2double(f(findToken(1)+1:findToken(2)-1));
    end
end
numList(isnan(numList)) = [];
numList = sort(numList);
numElements = length(numList);

Repeats = 1;

%initialise
% initFit_Collection_Plastic = NaN(Repeats * numElements,numDirs);
% fitDiff_Collection_Plastic = NaN(Repeats * numElements,4,numDirs);
initFit_Collection_Plastic = [];
fitDiff_Collection_Plastic = [];

delete(gcp('nocreate'));
myCluster = parcluster('local');
myCluster.NumWorkers = 10;
saveProfile(myCluster);
parpool(myCluster);

for i=1:numElements
    
    %show progress %
    disp([sprintf('%.2f',(round(10000*i/numElements)/100)) '%']);
    
    %load Pop%
    tempPop = load([folderDirectory '\replicate_' num2str(numList(i)) '_Pop.mat']);
    tempPop = tempPop.Pop;
    
    for k = 1:Repeats
        %show progress %
        %disp([sprintf('%.2f',(round(10000*k/Repeats)/100)) '%']);
        parfor j = 1:10
            %compute
            [tempA,tempB] = evalEvolvability_Inference_v2(tempPop, scenario);
            fitDiff_Collection_Plastic = [fitDiff_Collection_Plastic; tempA];
            initFit_Collection_Plastic = [initFit_Collection_Plastic; tempB];
        end
        %save
        if (exist([folderDirectory '\evolvabilityDS_v2.mat'],'dir'))
            disp('test');
           load([folderDirectory '\evolvabilityDS_v2.mat']);
        else
            fitDiff = []; initFit = [];
        end
        fitDiff = [fitDiff; fitDiff_Collection_Plastic]; %#ok<AGROW>
        initFit = [initFit; initFit_Collection_Plastic]; %#ok<AGROW>
        save([folderDirectory '\evolvabilityDS_v2.mat'],'fitDiff','initFit');
    end
end

end

