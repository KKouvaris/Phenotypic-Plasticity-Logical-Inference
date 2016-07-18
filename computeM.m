function [ uPhens, mPhens, recChange ] = computeM( Pop, scenario, mode, rand_a, rand_m)
%Evaluate phenotypic traits under random (single point) mutations

%constants
sigma_m = 0.5;
N = length(Pop);
S = size(Pop(1).Y,1);
phenSize = size(Pop(1).Z,1);
epsilon = phenSize;
alpha = 0.2; %probability of chaning the network architecture
tau = 20-1; %developmental time-steps

numY = numel(Pop(1).Y);
numZ = numel(Pop(1).Z);
numy = numel(Pop(1).y);

mNum = size(rand_a,1);

%initialise
uPhens = NaN(N*mNum,phenSize);
mPhens = NaN(N*mNum,phenSize);
recChange = zeros(N*mNum,1);

if (any(scenario == [0 1]) || mode == 4) %no input
    
    %get environmental cue
    [~, env_cue, ~] = drawNewTarget(scenario, epsilon);
    cue_signal = [0*env_cue; zeros(S-epsilon,1)];
    
    %get adult phenotypes from the original population
    phensBefore = zeros(N,phenSize);
    for i=1:N
        %develop
        P = ones(S,1);
        for t_dev=1:tau,
            P = (1+exp(-Pop(i).Y.*Pop(i).y*P + cue_signal)).^-1;
        end
        tempP = P;
        P = (1+exp(-Pop(i).Y.*Pop(i).y*P + cue_signal)).^-1;
        error = norm(P-tempP,2);
        
        %check stability
        if (error>0.01)
            phensBefore(i,:) = -Inf;
        else
            phensBefore(i,:) = Pop(i).Z.*Pop(i).z*P(epsilon+1:end,:);
        end
    end
    
    for i=1:N %foreach individual
        %i
        for j=1:mNum %foreach mutation
            
            %create new individual - reset
            Y = Pop(i).Y;
            y = Pop(i).y;
            Z = Pop(i).Z;
            z = Pop(i).z;
            
            %mutate pop
            a = rand();
            if (a<alpha) %change network structure
                r_idx=randi(numY + numZ);
                if (r_idx <= numY) %change gene interaction
                    Y(r_idx) = ~Y(r_idx);
                    recChange((i-1)*mNum+j,1) = 1;
                else %change gene-trait interaction
                    r_idx = r_idx - numY;
                    Z(r_idx) = ~Z(r_idx);
                end
            else %change weights
                
                r_idx = randi(numY + numZ);
                if (r_idx <= numy) %change gene interaction
                    y(r_idx) = y(r_idx) + normrnd(0,sigma_m);
                    recChange((i-1)*mNum+j,1) = 1;
                else %change gene-trait interaction
                    r_idx = r_idx - numy;
                    z(r_idx) = z(r_idx) + normrnd(0,sigma_m);
                end
            end
            
            %develop
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-Y.*y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-Y.*y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01)
                mPhens((i-1)*mNum+j,:) = -Inf;
            else
                mPhens((i-1)*mNum+j,:) = Z.*z*P(epsilon+1:end,:);
            end
            uPhens((i-1)*mNum+j,:) = phensBefore(i,:);
        end
    end
else
    for i=1:N %foreach individual
        for j=1:mNum %foreach mutation
            
            %get environmental cue
            [~, env_cue, mask] = drawNewTarget(scenario, epsilon);
            
            %compatible with 3bits phenotypes only
            switch mode
                case 2
                    mask = [1; 0; 1]; %AC
                case 3
                    mask = [1; 0; 0]; %A Only
            end
            
            cue_signal = [mask.*env_cue; zeros(S-epsilon,1)];
            
            %develop phensBefore
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-Pop(i).Y.*Pop(i).y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-Pop(i).Y.*Pop(i).y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error<0.01)
                uPhens((i-1)*mNum+j,:) = Pop(i).Z.*Pop(i).z*P(epsilon+1:end,:);
            end
            
            %create new individual - reset
            Y = Pop(i).Y;
            y = Pop(i).y;
            Z = Pop(i).Z;
            z = Pop(i).z;
            
            %mutate pop
            if (rand_a(j)<alpha) %change network structure
                r_idx = rand_m(j);
                if (r_idx <= numY) %change gene interaction
                    Y(r_idx) = ~Y(r_idx);
                    recChange((i-1)*mNum+j,1) = 1;
                else %change gene-trait interaction
                    r_idx = r_idx - numY;
                    Z(r_idx) = ~Z(r_idx);
                end
            else %change weights
                r_idx = rand_m(j);
                if (r_idx <= numy) %change gene interaction
                    y(r_idx) = y(r_idx) + normrnd(0,sigma_m);
                    recChange((i-1)*mNum+j,1) = 1;
                else %change gene-trait interaction
                    r_idx = r_idx - numy;
                    z(r_idx) = z(r_idx) + normrnd(0,sigma_m);
                end
            end
            
            %develop
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-Y.*y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-Y.*y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01), continue; end
            mPhens((i-1)*mNum+j,:) = Z.*z*P(epsilon+1:end,:);
        end
    end
end

end



