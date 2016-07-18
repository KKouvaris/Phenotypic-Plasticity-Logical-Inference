function [ fitDiff, initFit_Th ] = evalEvolvability_Inference_v2(initPop, scenario)
%Evaluate evolvability under directional selection as per Draghi's paper

%constants
s = 1;
sigma_m = 0.5;
N = length(initPop);
S = size(initPop(1).Y,1);
phenSize = size(initPop(1).Z,1);
epsilon = phenSize;
alpha = 0.2; %probability of chaning the network architecture
tau = 20-1; %developmental time-steps
offset = 2^-.5;
t_max = 1;

numY = numel(initPop(1).Y);
numZ = numel(initPop(1).Z);
numy = numel(initPop(1).y);
numz = numel(initPop(1).z);

min_cue = -0.5; %min cue value
max_cue =  0.5; %max cue value

%directions = [1 1 1; 1 1 0; 1 -1 0; 0 1 1; 0 1 -1; 1 0 1; 1 0 -1; 1 0 0; 0 1 0; 0 0 1];
switch scenario
    case {7,8}
        directions = [1 -1 0; 1 0 1; 0 -1 1];
    case {9,10,12}
        directions = [1 1 0; 0 1 -1; 1 0 -1];

end
numDirs = size(directions,1);
%normalise
for i = 1:numDirs
    directions(i,:) = directions(i,:)/norm(directions(i,:),2);
end

Pop(N) = struct();
newPop(N) = struct();
fitDiff = NaN(4,numDirs);
initFit_Th = NaN(1,numDirs);

for theta = 1:numDirs

    %initialise population
    for i=1:N
        Pop(i).Y = initPop(i).Y; %unipolar binary interaction matrix
        Pop(i).y = initPop(i).y; %signed contineous interaction matrix
        Pop(i).Z = initPop(i).Z; %unipolar binary interaction matrix
        Pop(i).z = initPop(i).z; %signed contineous interaction matrix
        Pop(i).fitness = initPop(i).fitness; %fitness
    end
    
    %get environmental cue
    env_cue = directions(theta,:);
    env_cue(env_cue>0) = max_cue;
    env_cue(env_cue<0) = min_cue;
    cue_signal = [env_cue zeros(1,S-epsilon)]';
    
    %evaluate mean initial fitness (t=0)
    fit_vector = zeros(N,1);
    for i=1:N
        %develop
        P = ones(S,1);
        for t_dev=1:tau,
            P = (1+exp(-Pop(i).Y.*Pop(i).y*P+cue_signal)).^-1;
        end
        tempP = P;
        P = (1+exp(-Pop(i).Y.*Pop(i).y*P+cue_signal)).^-1;
        error = norm(P-tempP,2);
        
        %check stability and calculate fitness for the respective angle
        if (error>0.01)
            Pop(i).fitness = NaN;
        else
            P = Pop(i).Z.*Pop(i).z*P(epsilon+1:end,:);
            d = sum((P-offset).*directions(theta,:)');
            Pop(i).fitness = (1+s)^d;
        end
        fit_vector(i,1) = Pop(i).fitness;
    end

    %remove unstable individuals from calculations
    fit_vector = fit_vector(~isnan(fit_vector));
    
    initFit = mean(fit_vector);
    initFit_Th(1,theta) = initFit;
    
    for t=1:t_max
        %evaluate weights
        fit_weights = zeros(N,1);
        for i=1:N, fit_weights(i,1) = Pop(i).fitness; end
        fit_weights(isnan(fit_weights)) = 0;
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
                    %accumulation(end)
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
                    P = (1+exp(-newPop(i).Y.*newPop(i).y*P+cue_signal)).^-1;
                end
                tempP = P;
                P = (1+exp(-newPop(i).Y.*newPop(i).y*P+cue_signal)).^-1;
                error = norm(P-tempP,2);
                
                %check stability
                if (error>0.01), continue; end
                
                P = newPop(i).Z.*newPop(i).z*P(epsilon+1:end,:);
                
                %fitness
                d = sum((P-offset).*directions(theta,:)');
                newPop(i).fitness = (1+s)^d;
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
        
        if (t==1 || t==10 || t==30 || t==100)
            if (t==1), ind_t = 1; end
            if (t==10), ind_t = 2; end
            if (t==30), ind_t = 3; end
            if (t==100), ind_t = 4; end
            fitDiff(ind_t,theta) = mean(([Pop(:).fitness]));
        end
    end
end
end

