function [ ] = evalAvgRN( scenario, mode )

%Scenario
% 6- {AB-} => AB-
% 7- {AB-,C neutral} => {AC0, BC0}
% 8- {AB-, BC-} => AC+
% 9- {AB+, BC-} => AC-
%10- {AB+, BC+} => AC+
%11- {AB+, AB-} => AC0
%12- {AB+, {BC+,BC-}=>BC0 } => AC0

%Modes
%1-for novel environment A,
%2-for novel environment B,
%3-for novel environment C

for iter=1:5

%Determine subFolder 
switch scenario
    case 6
        subFolder = ['Anticorrelation\' num2str(iter*50000)];      
    case 7
        subFolder = ['Zerocorrelation\' num2str(iter*50000)];              
    case 8
        subFolder = ['Inference_Cues_1\' num2str(iter*50000)];         
    case 9
        subFolder = ['Inference_Cues_2\' num2str(iter*50000)];      
    case 10
        subFolder = ['Inference_Cues_3\' num2str(iter*50000)];      
    case 11
        subFolder = ['Inference_Cues_4\' num2str(iter*50000)];          
    case 12
        subFolder = ['Inference_Cues_5\' num2str(iter*50000)];                   
end

N = 10;

%Define training set
min_phi =   (2*2^0.5)^-1;
max_phi = 3*(2*2^0.5)^-1;
min_cue = -.5;
max_cue =  .5;
Training_Sample_Size = 100;
Training_Set_Output = 0:1/(Training_Sample_Size-1):1;
switch scenario
    case 9
        Training_Set_Output = [ Training_Set_Output*(max_phi - min_phi) + min_phi;
            Training_Set_Output*(max_phi - min_phi) + min_phi;
            (1-Training_Set_Output)*(max_phi - min_phi) + min_phi];
end

Training_Set_Input = (Training_Set_Output' - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;

%Initialise
rNorm = [];

for i=10
    i
    load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\' subFolder '\replicate_' num2str(i) '_Pop.mat']);    
    rNorm = [rNorm evalReactionNorm(Pop,scenario)];
end

switch mode
    case 1
        focal_ = 1;
end

fig = figure;
for j=1:3
    subplot(1,3,j);
    if j ==2, title(['t = ' num2str(iter*50000)]); end
    hold on;
    plot(Training_Set_Input(:,focal_),Training_Set_Output(j,:)','r--','LineWidth',2);
    %for i=1:N, plot(Training_Set_Input(:,1),rNorm(:,(i-1)*3+j),'k--'); end
    plot(Training_Set_Input(:,1),rNorm(:,j),'k--'); 
    hold off;
    xlim([-.6 .6]);
    ylim([0 1.5]);
    axis square;
    xlabel('Env. Cue C');
    ylabel(['Trait ' num2str(j)]);
end

end

