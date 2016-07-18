function [] = evalMutCorr( scenario, mode, display, show, savefig, N, t_max, cores)

%Scenario
% 0- static environment,
% 1- heterogeneous environment,
% 2- plasticity,
% 3- plasticity binary,
% 4- plasticity test set,
% 5- static environment with env cue
% 6- {AB-} => AB-
% 7- {AB-,C neutral} => {AC0, BC0}
% 8- {AB-, BC-} => AC+
% 9- {AB+, BC-} => AC-
%10- {AB+, BC+} => AC+
%11- {AB+, AB-} => AC0
%12- {AB+, {BC+,BC-}=>BC0 } => AC0
%13- {ABC: A><B><C} => {AC0,AB0,BC0}
%14- {A,B,C} => {AC0,AB0,BC0}
%15- {ABC: AB+, BC-, AC-} static
%16-plastic as 9, non-linear FF, binary
%17-plastic as 9, non-linear FF, bipolar binary
%18-plastic, {ABCD: AB+ BC-}, S=20, [min_cue,max_cue]=[-1,1] change every generation

%Modes
%1-for past environments, 
%2-for novel environment AC,
%3-for novel environment A,
%4-no Input

%Display modes
%0-Evaluate Mutational Correlations
%1-Display Mutational Correlations

%Show modes
%0-Do not show figures
%1-Show figures

%Save modes
%0-Do not save figures
%1-Save figures

%cores
%number of used cores

%Arguments validation
if ~exist('t_max','var') || isempty(t_max)
    t_max = 50000; %default
end
if ~any( mode == 1:4) || ~any( display == [0 1]) || ~any( show == [0 1])
    warning('Error: Invalid arguments');
    return;
end

%Determine subFolder
folderDirectory = findFolder(scenario, N, t_max);
if ~exist(folderDirectory, 'dir')
    warning('Error: Directory does not exist');
    return; 
end

%Determine filename subString
switch mode
    case 2
        subString = '_AC';
    case 3
        subString = '_A';
    case 4
        subString = '_0';
    otherwise
        subString = '_default';
end

%get values
try
    load([folderDirectory '\replicate_1_Pop.mat']);
    numY = numel(Pop(1).Y);
    numZ = numel(Pop(1).Z);
    phenSize = size(Pop(1).Z,1); 
    corNum = nchoosek(phenSize,2); %number of pairwise correlations
catch
    warning(['Error: The following directory "' folderDirectory '\replicate_1_Pop.mat' '" does not exist']);
    return;
end

mNum = 300; %number of mutations for each individual

if (~display)
    %//////////// EVALUATE Unmutated and Mutated Phenotypes //////////%
    
    %get random numbers
    randDir = [folderDirectory '\rand_nums.mat'];
    if exist(randDir,'dir')
        %load
        load(randDir);
    else
        %create and save
        rand_a = rand(mNum,1);
        rand_m = randi(numY + numZ, mNum, 1);
        save([folderDirectory '\rand_nums.mat'],'rand_a','rand_m'); 
    end
    
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
    
    for i = 1:numElements
        
        %show progress %
        disp([sprintf('%.2f',(round(10000*i/numElements)/100)) '%']);
        
        %load Pop%
        load([folderDirectory '\replicate_' num2str(numList(i)) '_Pop.mat']);
        
        %evaluate unmutated and mutated phenotypes
        [uPhens,mPhens,recChange] = computeM(Pop,scenario,mode,rand_a,rand_m); %#ok<ASGLU>
        
        %save unmutated and mutated phenotypes
        save([folderDirectory '\replicate_' num2str(numList(i)) '_mutational_effects' subString '.mat'],'uPhens','mPhens','recChange');
        
    end

else
    %/////////////////////////DISPLAY M //////////////////////////////

    %determine the number of replicates in the respective scenario
    dirList = dir([folderDirectory '\*.mat']);
    numList = NaN(length(dirList),1);
    
    for id = 1:length(dirList)
        % Get the filename (minus the extension)
        [~, f] = fileparts(dirList(id).name);
        if ~isempty(strfind(f,['_mutational_effects' subString]))
            findToken = strfind(f,'_');
            numList(id) = str2double(f(findToken(1)+1:findToken(2)-1));
        end
    end
    
    numList(isnan(numList)) = [];
    numList = sort(numList);
    numElements = length(numList);
    
    %initialise
    corComb = nchoosek(1:phenSize,2); %get all possible pairwise combinations    
    meanCorr = NaN(numElements, corNum);
    meanVar = NaN(1, corNum);
    
    for i=1:numElements
        
        %show progress %
        disp([sprintf('%.2f',(round(10000*i/numElements)/100)) '%']);
        
        %load Pop%
        load([folderDirectory '\replicate_' num2str(numList(i)) '_mutational_effects' subString '.mat']);
        B = mPhens - uPhens;         %#ok<NODEF>
        
        %calculate mutational correlations
        
        %initialise
        tempCorr = NaN(N,corNum);
        
        for k=1:N %foreach individual
            tempB = B((k-1)*mNum+1:k*mNum,:);
            tempB = tempB(all(~isinf(tempB) & ~isnan(tempB),2),:); %filter
            if isempty(tempB), continue; end
            for cl=1:corNum
                tempComb = tempB(:,corComb(cl,:));
                tempCorr(k,cl) = corr(tempComb(:,1),tempComb(:,2));
            end
        end
        
        meanCorr(i,:) = mean(tempCorr(all(~isinf(tempCorr) & ~isnan(tempCorr),2),:));
        
        %calculate mutational variance
        tempB = B(all(~isinf(B) & ~isnan(B),2),:);
        meanVar(i) = mean(var(tempB));
        
        if (show || savefig)
            
            if (show), dispFig = 'on'; else dispFig = 'off'; end
            
            fig = figure('Visible',dispFig);

            for cl=1:corNum
                subplot(1,corNum,cl);
                hold on;
                title(['Corr:' num2str(meanCorr(i,cl))]);
                plot(B(logical(~recChange),corComb(cl,1)),B(logical(~recChange),corComb(cl,2)),'.'); %#ok<NODEF>
                plot(B(logical( recChange),corComb(cl,1)),B(logical( recChange),corComb(cl,2)),'.','color',[1 .5 0]);
                hold off;
                xlim([-1 1]); ylim([-1 1]); axis square;
                xlabel(['Trait ' num2alph(corComb(cl,1))]); 
                ylabel(['Trait ' num2alph(corComb(cl,2))]);              
            end
            
            if (savefig)
                print(fig,[folderDirectory '\mut_distribution_' subString '_' num2str(numList(i))],'-dpng');
            end
            
        end
        
    end
    meanCorr
    %Display results
    disp(['Var : ' num2str(mean(meanVar))]);
    for cl=1:corNum
       disp(['Corr (' num2alph(corComb(cl,1)) ',' num2alph(corComb(cl,2)) '): ' num2str(mean(meanCorr(:,cl)))]);
    end   
end