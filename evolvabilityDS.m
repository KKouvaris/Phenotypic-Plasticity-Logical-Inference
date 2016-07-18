function [ initFit_Collection_Plastic, fitDiff_Collection_Plastic ] = evolvabilityDS( scenario, t_max )

folderDirectory = findFolder(scenario, t_max);

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

Repeats = 3;
numDirs = 6;

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
    
    for k = 1:Repeats
        %show progress %
        disp([sprintf('%.2f',(round(10000*k/Repeats)/100)) '%']);
        parfor j = 1:10
            %compute
            [tempA,tempB] = evalEvolvability_Inference(tempPop.Pop, scenario);
            fitDiff_Collection_Plastic = [fitDiff_Collection_Plastic; tempA];
            initFit_Collection_Plastic = [initFit_Collection_Plastic; tempB];
            %[fitDiff_Collection_Plastic((i-1)*Repeats+j,:,:),initFit_Collection_Plastic((i-1)*Repeats+j,:)] = evalEvolvability_Inference(Pop, scenario);
        end
        %save
        if (exist([folderDirectory '\evolvabilityDS.mat'],'dir'))
            load([folderDirectory '\evolvabilityDS.mat']);
        else
            fitDiff = []; initFit = [];
        end
        fitDiff = [fitDiff; fitDiff_Collection_Plastic]; %#ok<AGROW>
        initFit = [initFit; initFit_Collection_Plastic]; %#ok<AGROW>
        save([folderDirectory '\evolvabilityDS.mat'],'fitDiff','initFit');
    end
end

% plastic_vector = NaN(4,numDirs);
%
% for i=1:4
%     temp_plastic = mean(log(fitDiff_Collection_Plastic(:,i,:)),1);
%
%     for j=1:numDirs
%         plastic_vector(i,j) = temp_plastic(:,:,j);
%     end
% end
%
% fig = figure;
% for i=1:4
%     plastic_vector(i,:)
%     mean(initFit_Collection_Plastic,1)
%     plastic_vector(i,:) = plastic_vector(i,:) - mean(log(initFit_Collection_Plastic),1);
%
%     subplot(2,2,i);
%     if i==1, title('t = 1'); end
%     if i==2, title('t = 10'); end
%     if i==3, title('t = 30'); end
%     if i==4, title('t = 100'); end
%     hold on; plot(plastic_vector(i,:),'ko'); xlabel('Directions'); ylabel('log-fitness difference'); xlim([0 numDirs+1]);
% end

end

