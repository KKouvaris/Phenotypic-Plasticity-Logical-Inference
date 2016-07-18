function [ M ] = computeM( Pop, scenario )
%COMPUTEM Summary of this function goes here
%   Detailed explanation goes here

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
numz = numel(Pop(1).z);

min_phi =   (2*2^0.5)^-1;
max_phi = 3*(2*2^0.5)^-1;

min_cue = -0.5;
max_cue =  0.5;

mNum = 100; %number of total mutations, sample size

if (scenario == 2)
	phi = rand()*(max_phi - min_phi) + min_phi;
	env_cue = (phi - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;
	cue_signal = [env_cue*ones(1,epsilon) zeros(1,S-epsilon)]';    
else
    env_cue = 0;
    cue_signal = [env_cue*ones(1,epsilon) zeros(1,S-epsilon)]';
end    
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
            continue;
        end
        phensBefore(i,:) = Pop(i).Z.*Pop(i).z*P(epsilon+1:end,:);
    end
    
    M = zeros(phenSize);
    M_Variance = zeros(1,phenSize);
    M_Cov = zeros(1,phenSize);
    phenDiff_Sum = zeros(N,phenSize);
    for j=1:mNum %foreach mutation
        j
        %create new pop - mutants
        newPop = struct();
        for i=1:N
            newPop(i).Y = Pop(i).Y;
            newPop(i).y = Pop(i).y;
            newPop(i).Z = Pop(i).Z;
            newPop(i).z = Pop(i).z;
        end
        
        %mutate pop
        if (rand()<alpha) %change network structure
            r_idx = randi(numY + numZ);
            if (r_idx <= numY) %change gene interaction
                for i=1:N, newPop(i).Y(r_idx) = ~newPop(i).Y(r_idx); end
            else %change gene-trait interaction
                r_idx = r_idx - numY;
                for i=1:N, newPop(i).Z(r_idx) = ~newPop(i).Z(r_idx); end
            end
        else %change weights
            r_idx = randi(numy + numz);
            if (r_idx <= numy) %change gene interaction
                for i=1:N, newPop(i).y(r_idx) = newPop(i).y(r_idx) + normrnd(0,sigma_m); end
            else %change gene-trait interaction
                r_idx = r_idx - numy;
                for i=1:N, newPop(i).z(r_idx) = newPop(i).z(r_idx) + normrnd(0,sigma_m); end
            end
        end

        phensAfter = zeros(N,phenSize);
        %develop
        for i=1:N
            P = ones(S,1);
            for t_dev=1:tau,
                P = (1+exp(-newPop(i).Y.*newPop(i).y*P + cue_signal)).^-1;
            end
            tempP = P;
            P = (1+exp(-newPop(i).Y.*newPop(i).y*P + cue_signal)).^-1;
            error = norm(P-tempP,2);
            
            %check stability
            if (error>0.01)
                phensAfter(i,:) = -Inf;
                continue; 
            end
            phensAfter(i,:) = newPop(i).Z.*newPop(i).z*P(epsilon+1:end,:);
        end
        
        phenDiff = phensAfter-phensBefore;
        %phenDiff_Sum = phenDiff_Sum + phenDiff;
        M_Variance = M_Variance + sum((phenDiff).^2)/N;
        M_Cov = M_Cov + sum(phensBefore(:,1).*phensBefore(:,2)-phensAfter(:,1).*phensAfter(:,2))/N;
        sum(phensBefore(:,1).*phensBefore(:,2)-phensAfter(:,1).*phensAfter(:,2))/N
        %M = M + phenDiff'*phenDiff/N; %cov(phenDiff(:,1),phenDiff(:,2));
        M = M + cov(phensAfter-phensBefore);
    end
    M = M
    M_Variance = M_Variance/mNum
    M_Cov = M_Cov/mNum
    %phenDiff_Sum


end

