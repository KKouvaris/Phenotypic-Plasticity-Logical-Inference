function [ ] = plotEvolDS( scenario_A, scenario_B, t_max_A, t_max_B )

folderDirectory_A = findFolder(scenario_A, t_max_A);
folderDirectory_B = findFolder(scenario_B, t_max_B);

%load fitDiff and initFit files
load([folderDirectory_A '\evolvabilityDS.mat']);
fitDiff_A = fitDiff;
initFit_A = initFit;    
load([folderDirectory_B '\evolvabilityDS.mat']);
fitDiff_B = fitDiff;
initFit_B = initFit; 

numDirs = size(initFit,2);

fig = figure;
for i=1:4
    
    temp_diffs_A = log(fitDiff_A(i:4:end,:))-log(initFit_A);
    temp_diffs_B = log(fitDiff_B(i:4:end,:))-log(initFit_B);
    
    size(temp_diffs_A)

    temp_diffs_A = (fitDiff_A(i:4:end,:))-(initFit_A);
    temp_diffs_B = (fitDiff_B(i:4:end,:))-(initFit_B);    
    vector_A = temp_diffs_A;
    vector_B = temp_diffs_B;    
    
%     vector_A = [];
%     vector_B = [];
%     for j=1:20
%         vector_A = [vector_A; mean(temp_diffs_A((j-1)*30+1:j*30,:),1)];
%         vector_B = [vector_B; mean(temp_diffs_B((j-1)*30+1:j*30,:),1)];
%     end
    
%     mean_vector = mean(temp_diffs,1);
%     std_vector = std(temp_diffs);
    


    subplot(2,2,i);
    switch i
        case 1
            title('t = 1')
        case 2
            title('t = 10');
        case 3
            title('t = 30');
        case 4
            title('t = 100');
    end
    
    hold on; plot(mean(vector_A),'rx'); plot(mean(vector_B),'bo'); xlabel('Directions'); ylabel('log-fitness difference'); xlim([0 numDirs+1]);
end
end

