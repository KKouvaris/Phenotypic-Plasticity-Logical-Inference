function [rNorm] = evalReactionNorm( Pop, scenario )

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

%constants
sigma_m = 0.5;
N = length(Pop);
S = size(Pop(1).Y,1);
phenSize = size(Pop(1).Z,1);
epsilon = phenSize;
omega = 0.2; %selection pressure
alpha = 0.2; %probability of chaning the network architecture
tau = 20-1; %developmental time-steps

min_phi =   (2*2^0.5)^-1;
max_phi = 3*(2*2^0.5)^-1;

min_cue = -.5;
max_cue =  .5;

% min_cue = -1;
% max_cue =  1;

% a = (max_phi - min_phi) * (min_cue+0.5) + min_phi;
% b = (max_phi - min_phi) * (max_cue+0.5) + min_phi;
% 
% min_phi = a;
% max_phi = b;

%Define training set
Training_Sample_Size = 100;
Training_Set_Output = 0:1/(Training_Sample_Size-1):1;
switch scenario
    case {9,16}
        Training_Set_Output = [ Training_Set_Output*(max_phi - min_phi) + min_phi;
            Training_Set_Output*(max_phi - min_phi) + min_phi;
            (1-Training_Set_Output)*(max_phi - min_phi) + min_phi];
    case 13
        Training_Set_Output = [ Training_Set_Output*(max_phi - min_phi) + min_phi;
            Training_Set_Output*(max_phi - min_phi) + min_phi;
            Training_Set_Output*(max_phi - min_phi) + min_phi];      
    case 17
        Training_Set_Output = [ 2*Training_Set_Output-1;
            2*Training_Set_Output-1;
            2*(1-Training_Set_Output)-1];        
end

Training_Set_Input = (Training_Set_Output' - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;

if (scenario == 17), Training_Set_Input = Training_Set_Output'; end

%initialise
rNorm = zeros(Training_Sample_Size,phenSize);

%calculate maximal fitness over the training set
training_fit = zeros(N,1);
for ind = 1:N
    sum_fit = 0;
    for sample = 1:Training_Sample_Size
        
        cue_signal = [ Training_Set_Input(sample,:)  zeros(1,S-epsilon)]';
        
        %develop
        P = ones(S,1);
        for t_dev=1:tau,
            P = (1+exp(-Pop(ind).Y.*Pop(ind).y*P + cue_signal)).^-1;
        end
        tempP = P;
        P = (1+exp(-Pop(ind).Y.*Pop(ind).y*P + cue_signal)).^-1;
        error = norm(P-tempP,2);
        
        %check stability and evaluate fitness
        if (error<0.01)
            P = Pop(ind).Z.*Pop(ind).z*P(epsilon+1:end,:);
            if (scenario == 16), P = tanh(P); end
            d = norm(P(:)-Training_Set_Output(:,sample),2);
            sum_fit = sum_fit + exp(-(d^2)/(2*omega));
        end
    end
    training_fit(ind,1) = sum_fit/Training_Sample_Size;
    
end
[~, max_index] = max(training_fit)


%max_index = 652;

%calculate reaction norm wtr A
for sample = 1:Training_Sample_Size
    
    cue_signal = [ Training_Set_Input(sample,:) zeros(1,S-epsilon)]';
    
    %develop
    P = ones(S,1);
    for t_dev=1:tau,
        P = (1+exp(-Pop(max_index).Y.*Pop(max_index).y*P + cue_signal)).^-1;
    end
    tempP = P;
    P = (1+exp(-Pop(max_index).Y.*Pop(max_index).y*P + cue_signal)).^-1;
    error = norm(P-tempP,2);
    
    %check stability and evaluate fitness
    if (error<0.01)
        P = Pop(max_index).Z.*Pop(max_index).z*P(epsilon+1:end,:);
        if (scenario == 16), P = tanh(P); end
        rNorm(sample,:) = P;
    end
end

fig = figure;
for j=1:3
    subplot(1,3,j);
    hold on;
    plot(Training_Set_Input(:,1),Training_Set_Output(j,:)','r--','LineWidth',2);
    plot(Training_Set_Input(:,1),rNorm(:,j),'k--');
    hold off;
    %xlim([-1.1 1.1]);
    %ylim([0 1.5]);
    axis square;
    xlabel('Environmental Cues');
    ylabel(['Trait ' num2str(j)]);
end

%print(fig,['reaction norm' num2str(time_label)],'-dpng');

end
