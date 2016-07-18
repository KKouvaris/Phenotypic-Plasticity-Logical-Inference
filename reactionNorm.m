N = 1000; %population size;

modules = 1; %number of modules
S = modules * 20; %number of genes
epsilon = modules * 2; %environmental cues
phenSize = epsilon; %number of phenotypic traits
alpha = 0.2; %probability of chaning the network architecture
sigma_m = 0.5; %std variation
omega = 0.2; %selection pressure
tau = 20-1; %developmental time-steps

min_phi =   (2*2^0.5)^-1;
max_phi = 3*(2*2^0.5)^-1;

min_cue = -1;
max_cue =  1;

a = (max_phi - min_phi) *(min_cue+0.5)+min_phi;
b = (max_phi - min_phi) *(max_cue+0.5)+min_phi;

min_phi = a;
max_phi = b;

Training_Sample_Size = 100;
Training_Set_Output = [0:1/(Training_Sample_Size-1):1]*(max_phi - min_phi) + min_phi;
Training_Set_Input = (Training_Set_Output' - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;

for time=[1 10 20 50 51 52]
    
    %time=52
    
    iter = 1;
    rNorm_1 = zeros(Training_Sample_Size,iter);
    rNorm_2 = zeros(Training_Sample_Size,iter);
    
    for i=1:iter
        
        i
        
        %%load Pop%
        %load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50k\2outputs\replicate_' num2str(i) '_Pop_Array.mat']);
        
        Pop = Pop_Array{time};
        
        
        %find max fitness
        %calculate fitness over the training set
%         training_fit = zeros(N,1);
%         for ind = 1:N
%             sum_fit = 0;
%             for sample = 26:75 %1:Training_Sample_Size
%                 
%                 cue_signal = [Training_Set_Input(sample,1)*ones(1,epsilon) zeros(1,S-epsilon)]';
%                 
%                 %develop
%                 P = ones(S,1);
%                 for t_dev=1:tau,
%                     P = (1+exp(-Pop(ind).Y.*Pop(ind).y*P + cue_signal)).^-1;
%                 end
%                 tempP = P;
%                 P = (1+exp(-Pop(ind).Y.*Pop(ind).y*P + cue_signal)).^-1;
%                 error = norm(P-tempP,2);
%                 
%                 %check stability and evaluate fitness
%                 if (error<0.01)
%                     P = Pop(ind).Z.*Pop(ind).z*P(epsilon+1:end,:);
%                     d = norm(P-Training_Set_Output(1,sample)*ones(size(P)),2);
%                     sum_fit = sum_fit + exp(-(d^2)/(2*omega));
%                 end
%             end
%             training_fit(ind,1) = sum_fit/Training_Sample_Size;
%             
%         end
%         [~, max_index] = max(training_fit);
        max_index = 1;
        
        %calculate reaction norm
        for sample = 1:Training_Sample_Size
            
            cue_signal = [Training_Set_Input(sample,1)*ones(1,epsilon) zeros(1,S-epsilon)]';
            
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
                rNorm_1(sample,i) = P(1);
                rNorm_2(sample,i) = P(2);
            end
        end
    end
    
    time_label = time;
    if (time == 51), time_label = 10000; end
    if (time == 52), time_label = 50000; end
    
    fig = figure;
    subplot(1,2,1);
    title(['Generation ' num2str(time_label)]);
    hold on;
    plot(Training_Set_Input,Training_Set_Output,'r--','LineWidth',2);
    for i=1:iter, plot(Training_Set_Input,rNorm_1(:,i),'k--'); end;
    hold off;
    %xlim([-1.1 1.1]);
    %ylim([-1.5 1.5]);
    axis square;
    xlabel('Environmental Cues');
    ylabel('Trait 1');
    subplot(1,2,2);
    title(['Generation ' num2str(time_label)]);
    hold on;
    plot(Training_Set_Input,Training_Set_Output,'r--','LineWidth',2);
    for i=1:iter, plot(Training_Set_Input,rNorm_2(:,i),'k--'); end;
    hold off;
    %xlim([-1.1 1.1]);
    %ylim([-1.5 1.5]);
    axis square;
    xlabel('Environmental Cues');
    ylabel('Trait 2');
    
    print(fig,['reaction norm' num2str(time_label)],'-dpng');
end
