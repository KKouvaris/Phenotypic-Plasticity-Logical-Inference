function [ Parents_Traits, Child_Traits ] = computeG( Pop, scenario )
%COMPUTEM Summary of this function goes here
%   Detailed explanation goes here

%constants
N = length(Pop);
S = size(Pop(1).Y,1);
phenSize = size(Pop(1).Z,1);
epsilon = phenSize;
tau = 20-1; %developmental time-steps

numY = numel(Pop(1).Y);
numZ = numel(Pop(1).Z);

min_phi =   (2*2^0.5)^-1;
max_phi = 3*(2*2^0.5)^-1;

min_cue = -0.5;
max_cue =  0.5;

nEnv = 10;

if (scenario == 2)
    
    C = combnk(1:N,2); %get parents combinations
    Parents_Traits = NaN(size(C,1)*nEnv,2);
    Child_Traits = NaN(size(C,1)*nEnv,2);
    
    for env = 1:nEnv
        
        %choose env
        phi = rand()*(max_phi - min_phi) + min_phi;
        phi = kron(phi,ones(2,1));
        env_cue = (phi - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;
        cue_signal = [env_cue; zeros(S-epsilon,1)];
        
        env
        
        for i=1:size(C,1)
            %[i size(C,1)]
            idx_parent_1 = C(i,1);
            idx_parent_2 = C(i,2);
            
            %create child
            Y = Pop(idx_parent_1).Y;
            y = Pop(idx_parent_1).y;
            Z = Pop(idx_parent_1).Z;
            z = Pop(idx_parent_1).z;
            
            %crossover
            idx_rnd_Y_rows = logical(randi(2,S,1)-1);
            Y(idx_rnd_Y_rows,:) = Pop(idx_parent_2).Y(idx_rnd_Y_rows,:);
            y(idx_rnd_Y_rows,:) = Pop(idx_parent_2).y(idx_rnd_Y_rows,:);
            
            idx_rnd_Z_elem = logical(randi(2,2,S-epsilon)-1);
            Z(idx_rnd_Z_elem) = Pop(idx_parent_2).Z(idx_rnd_Z_elem);
            z(idx_rnd_Z_elem) = Pop(idx_parent_2).z(idx_rnd_Z_elem);
            
            %develop
            
            %parent 1
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-Pop(idx_parent_1).Y.*Pop(idx_parent_1).y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-Pop(idx_parent_1).Y.*Pop(idx_parent_1).y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01), continue; end
            
            P_Parent_1 = Pop(idx_parent_1).Z.*Pop(idx_parent_1).z*P(epsilon+1:end,:);
            
            %parent 2
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-Pop(idx_parent_2).Y.*Pop(idx_parent_2).y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-Pop(idx_parent_2).Y.*Pop(idx_parent_2).y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01), continue; end
            
            P_Parent_2 = Pop(idx_parent_2).Z.*Pop(idx_parent_2).z*P(epsilon+1:end,:);
            
            %child
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-Y.*y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-Y.*y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01), continue; end
            
            P_Child = Z.*z*P(epsilon+1:end,:);
            P_Parent = mean([P_Parent_1 P_Parent_2],2);
            
            Parents_Traits((env-1)*nEnv+i,:) = P_Parent;
            Child_Traits((env-1)*nEnv+i,:) = P_Child;
        end
    end
else
    C = combnk(1:N,2); %get parents combinations
    Parents_Traits = NaN(size(C,1),2);
    Child_Traits = NaN(size(C,1),2);
    for i=1:size(C,1)
        %[i size(C,1)]
        idx_parent_1 = C(i,1);
        idx_parent_2 = C(i,2);
        
        %create child
        Y = Pop(idx_parent_1).Y;
        y = Pop(idx_parent_1).y;
        Z = Pop(idx_parent_1).Z;
        z = Pop(idx_parent_1).z;
        
        %crossover
        idx_rnd_Y_rows = logical(randi(2,S,1)-1);
        Y(idx_rnd_Y_rows,:) = Pop(idx_parent_2).Y(idx_rnd_Y_rows,:);
        y(idx_rnd_Y_rows,:) = Pop(idx_parent_2).y(idx_rnd_Y_rows,:);
        
        idx_rnd_Z_elem = logical(randi(2,2,S-epsilon)-1);
        Z(idx_rnd_Z_elem) = Pop(idx_parent_2).Z(idx_rnd_Z_elem);
        z(idx_rnd_Z_elem) = Pop(idx_parent_2).z(idx_rnd_Z_elem);
        
        %develop
        
        %parent 1
        P = ones(S,1);
        for t_dev=1:tau,
            P = (1+exp(-Pop(idx_parent_1).Y.*Pop(idx_parent_1).y*P)).^-1;
        end
        tempP = P;
        P = (1+exp(-Pop(idx_parent_1).Y.*Pop(idx_parent_1).y*P)).^-1;
        error = norm(P-tempP,2);
        
        %check stability
        if (error>0.01), continue; end
        
        P_Parent_1 = Pop(idx_parent_1).Z.*Pop(idx_parent_1).z*P(epsilon+1:end,:);
        
        %parent 2
        P = ones(S,1);
        for t_dev=1:tau,
            P = (1+exp(-Pop(idx_parent_2).Y.*Pop(idx_parent_2).y*P)).^-1;
        end
        tempP = P;
        P = (1+exp(-Pop(idx_parent_2).Y.*Pop(idx_parent_2).y*P)).^-1;
        error = norm(P-tempP,2);
        
        %check stability
        if (error>0.01), continue; end
        
        P_Parent_2 = Pop(idx_parent_2).Z.*Pop(idx_parent_2).z*P(epsilon+1:end,:);
        
        %child
        P = ones(S,1);
        for t_dev=1:tau,
            P = (1+exp(-Y.*y*P)).^-1;
        end
        tempP = P;
        P = (1+exp(-Y.*y*P)).^-1;
        error = norm(P-tempP,2);
        
        %check stability
        if (error>0.01), continue; end
        
        P_Child = Z.*z*P(epsilon+1:end,:);
        P_Parent = mean([P_Parent_1 P_Parent_2],2);
        
        Parents_Traits(i,:) = P_Parent;
        Child_Traits(i,:) = P_Child;
    end
end

end


