function [saveObj] = GRN_Fun(scenario, init_t, init_Pop)

%Arguments validation

if ~exist('init_t','var') || isempty(init_t) || mod(init_t,50000) ~= 0
    init_t = 0; %default
end
if init_t ~= 0 && (~exist('init_Pop','var') || isempty(init_Pop))
    warning('Error with initial population');
    return;
end
if ~exist('scenario','var') || isempty(scenario) || (~any( scenario == 0:18))
    warning('Error with scenario number');
    return;
end

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
%16-plastic as 9, non-linear FF, binary
%17-plastic as 9, non-linear FF, bipolar binary
%18-plastic, {ABCD: AB+ BC-}, S=20, [min_cue,max_cue]=[-1,1] change every generation

if (scenario == 17)
    N = 500; %population size;
else
    N = 1000; %population size;
end

modules = 1; %number of modules
if (scenario == 17)
    S = modules * 10; %number of genes
else
   S = modules * 20; %number of genes
end

epsilon = modules * 2; %default number of environmental cues
if (any( scenario == [7:10 12:17])), epsilon = 3; end
if (any( scenario == 18)), epsilon = 4; end
phenSize = epsilon; %number of phenotypic traits
alpha = 0.2; %probability of chaning the network architecture
sigma_m = 0.5; %std variation
omega = 0.2; %selection pressure
tau = 20-1; %developmental time-steps

t_max = 50000;

%Pop_Array = cell(52,1);

% Training_Sample_Size = 20;
% Training_Set_Output = [0:1/Training_Sample_Size:1]*(max_phi - min_phi) + min_phi;
% Training_Set_Input = (Training_Set_Output' - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;
%
% Test_Sample_Size = 20;
% Test_Set_Output = [0:1/Test_Sample_Size:1]*(max_phi - min_phi) + min_phi;
% Test_Set_Output(1:Test_Sample_Size/2) = Test_Set_Output(1:Test_Sample_Size/2) - (max_phi - min_phi)/2;
% Test_Set_Output(Test_Sample_Size/2+1:end) = Test_Set_Output(Test_Sample_Size/2+1:end) + (max_phi - min_phi)/2;
% Test_Set_Input = 2*((Test_Set_Output' - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue);

if (init_t == 0)
    
    %Initialise
    
    %find an initial GRN that produces a stable expression and fitness>0.00005
    found = 0;
    while(~found)
        %randomise individual
        founder_Y = randi(2,S)-1;
        founder_y = normrnd(0,sigma_m,S);
        founder_Z = randi(2,phenSize,S-epsilon)-1;
        founder_z = normrnd(0,sigma_m,phenSize,S-epsilon);
        
        %choose phenotypic target and environmental cue at random
        [phi, env_cue, mask] = drawNewTarget(scenario, epsilon);
        cue_signal = [mask.*env_cue; zeros(S-epsilon,1)];
        
        %develop
        P = ones(S,1);
        for i=1:tau,
            P = (1+exp(-founder_Y.*founder_y*P + cue_signal)).^-1;
        end
        tempP = P;
        P = (1+exp(-founder_Y.*founder_y*P + cue_signal)).^-1;
        error = norm(P-tempP,2);
        
        %check stability
        if (error<0.01)
            P = founder_Z.*founder_z*P(epsilon+1:end,:);
            if (any(scenario == 16:17)), P = tanh(P); end
        else
            continue;
        end
        
        %check fitness
        d = norm(mask.*(P-phi),2);
        init_fit = exp(-(d^2)/(2*omega));
        if init_fit > 0.00005, found = 1; end
    end
    %
    % %create a population of clones
    Pop(N) = struct();
    for i=1:N
        Pop(i).Y = founder_Y; %unipolar binary interaction matrix
        Pop(i).y = founder_y; %signed contineous interaction matrix
        Pop(i).Z = founder_Z; %unipolar binary interaction matrix
        Pop(i).z = founder_z; %signed contineous interaction matrix
        Pop(i).fitness = init_fit; %fitness
    end
    
else
    Pop = init_Pop;
end

newPop(N) = struct();
for i=1:N
    newPop(i).Y = Pop(i).Y; %unipolar binary interaction matrix
    newPop(i).y = Pop(i).y; %signed contineous interaction matrix
    newPop(i).Z = Pop(i).Z; %unipolar binary interaction matrix
    newPop(i).z = Pop(i).z; %signed contineous interaction matrix
    newPop(i).fitness = Pop(i).fitness; %fitness
end

numY = numel(Pop(i).Y);
numy = numel(Pop(i).y);
numZ = numel(Pop(i).Z);
numz = numel(Pop(i).z);

%output_b = zeros(t_max+1, 2*(S^2+2*(S-epsilon)));
%output_b(1,:) = [Pop(1).Y(:)' Pop(1).y(:)' Pop(1).Z(:)' Pop(1).z(:)'];
fit_output = zeros(t_max,3);
%     fit_training = zeros(t_max,3);
%     fit_test = zeros(t_max,3);

for t = 1:t_max
    
    %alternate selective environments
    if (~any(scenario == [0 5 15])) %non-static
        [phi, env_cue, mask] = drawNewTarget(scenario, epsilon);
        cue_signal = [mask.*env_cue; zeros(S-epsilon,1)];
    end
    
    %evaluate weights
    fit_weights = zeros(N,1);
    for i=1:N, fit_weights(i,1) = Pop(i).fitness; end
    fit_weights = fit_weights/mean(fit_weights);
    accumulation = cumsum(fit_weights);
    
    %create new population
    for i=1:N
        childFound = 0;
        while(~childFound)
            
            %select parents - roulette wheel selection
            %parent 1
            p = rand() * accumulation(end);
            idx_parent_1 = -1;
            for index = 1 : length(accumulation)
                if (accumulation(index) > p)
                    idx_parent_1 = index;
                    break;
                end
            end
            
            %parent 2
            while (true)
                p = rand() * accumulation(end);
                idx_parent_2 = -1;
                for index = 1 : length(accumulation)
                    if (accumulation(index) > p)
                        idx_parent_2 = index;
                        break;
                    end
                end
                if (idx_parent_1 ~= idx_parent_2), break; end
            end
            
            %create child
            newPop(i).Y = Pop(idx_parent_1).Y;
            newPop(i).y = Pop(idx_parent_1).y;
            newPop(i).Z = Pop(idx_parent_1).Z;
            newPop(i).z = Pop(idx_parent_1).z;
            newPop(i).fitness = -1;
            
            %crossover
            idx_rnd_Y_rows = logical(randi(2,S,1)-1);
            newPop(i).Y(idx_rnd_Y_rows,:) = Pop(idx_parent_2).Y(idx_rnd_Y_rows,:);
            newPop(i).y(idx_rnd_Y_rows,:) = Pop(idx_parent_2).y(idx_rnd_Y_rows,:);
            
            idx_rnd_Z_elem = logical(randi(2,2,S-epsilon)-1);
            newPop(i).Z(idx_rnd_Z_elem) = Pop(idx_parent_2).Z(idx_rnd_Z_elem);
            newPop(i).z(idx_rnd_Z_elem) = Pop(idx_parent_2).z(idx_rnd_Z_elem);
            
            %mutate pop
            mNum = poissrnd(0.02);
            for j=1:mNum %foreach mutation
                if (rand()<alpha) %change network structure
                    r_idx = randi(numY+numZ);
                    if (r_idx <= numY) %change gene interaction
                        newPop(i).Y(r_idx) = ~newPop(i).Y(r_idx);
                    else %change gene-trait interaction
                        r_idx = r_idx - numY;
                        newPop(i).Z(r_idx) = ~newPop(i).Z(r_idx);
                    end
                else %change weights
                    r_idx = randi(numy+numz);
                    if (r_idx <= numy) %change gene interaction
                        newPop(i).y(r_idx) = newPop(i).y(r_idx) + normrnd(0,sigma_m);
                    else %change gene-trait interaction
                        r_idx = r_idx - numy;
                        newPop(i).z(r_idx) = newPop(i).z(r_idx) + normrnd(0,sigma_m);
                    end
                end
            end
            
            %develop
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-newPop(i).Y.*newPop(i).y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-newPop(i).Y.*newPop(i).y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01), continue; end
            
            P = newPop(i).Z.*newPop(i).z*P(epsilon+1:end,:);
            if (any(scenario == 16:17)), P = tanh(P); end
            
            %fitness
            d = norm(mask.*(P-phi),2);
            newPop(i).fitness = exp(-(d^2)/(2*omega));
            childFound = 1;
        end
    end
    
    %copy new generation
    for i=1:N
        Pop(i).Y = newPop(i).Y;
        Pop(i).y = newPop(i).y;
        Pop(i).Z = newPop(i).Z;
        Pop(i).z = newPop(i).z;
        Pop(i).fitness = newPop(i).fitness;
    end
    
    %         if (t<=50), Pop_Array{t} = Pop; end
    %         if (t==10000), Pop_Array{51} = Pop; end
    %         if (t==t_max), Pop_Array{52} = Pop; end
    
    if (mod(t,50) == 0)
        
        %store fitness (min, mean, max)
        fit_vector = zeros(N,1);
        for ind = 1:N
            fit_vector(ind,1) = Pop(ind).fitness;
        end
        fit_output(t,:) = [min(fit_vector); mean(fit_vector); max(fit_vector)];
        
        %         %calculate fitness over the training set
        %         training_fit = zeros(N,1);
        %         for ind = 1:N
        %             sum_fit = 0;
        %             for sample = 1:Training_Sample_Size
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
        %         fit_training(t,:) = [min(training_fit); mean(training_fit); max(training_fit)];
        %
        %         %calculate fitness over the test set
        %         test_fit = zeros(N,1);
        %         for ind = 1:N
        %             sum_fit = 0;
        %             for sample = 1:Test_Sample_Size
        %
        %                 cue_signal = [Test_Set_Input(sample,1)*ones(1,epsilon) zeros(1,S-epsilon)]';
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
        %                     d = norm(P-Test_Set_Output(1,sample)*ones(size(P)),2);
        %                     sum_fit = sum_fit + exp(-(d^2)/(2*omega));
        %                 end
        %             end
        %             test_fit(ind,1) = sum_fit/Test_Sample_Size;
        %
        %         end
        %         fit_test(t,:) = [min(test_fit); mean(test_fit); max(test_fit)];
        
        disp([num2str(100*(t/t_max)) '%']);
        disp([min(fit_vector) mean(fit_vector) max(fit_vector)]);
        %         disp([min(training_fit) mean(training_fit) max(training_fit)]);
        %         disp([min(test_fit) mean(test_fit) max(test_fit)]);
    end
end

%save
saveObj = struct();
saveObj.Pop = Pop;
saveObj.t_max = t_max+init_t;
saveObj.N = N;

%saveMyVariables(scenario,modules*2,init_t + t_max,N,fit_output,Pop);


