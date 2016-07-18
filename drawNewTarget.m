function [ target, cue, mask ] = drawNewTarget( scenario, epsilon )
%choose new phenotypic target and environmental cue at random

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

%constants
min_phi =   (2*2^0.5)^-1; %min target value
max_phi = 3*(2*2^0.5)^-1; %max target value
min_cue = -0.5; %min cue value
max_cue =  0.5; %max cue value

if scenario == 18
    min_cue = -1; %min cue value
    max_cue =  1; %max cue value
end

mask_1 = [1; 1; 0];
mask_2 = [0; 1; 1];
mask_3 = [1; 1; 0; 0];
mask_4 = [0; 1; 1; 0];


%get new target
rnd = rand();
switch scenario
    case {0, 1, 2, 5}
        phi = [rnd; rnd];
    case 3
        rnd = round(rnd);
        phi = [rnd; rnd];
    case 6
        phi = [rnd; 1-rnd];
    case {7, 8}
        phi = [rnd; (1 - rnd); rnd];
    case {9, 15}
        phi = [rnd; rnd; (1 - rnd)];
    case 10
        phi = [rnd; rnd; rnd];
    case 11
        if (rand()<.5), phi = [rnd; rnd]; else phi = [rnd; (1 - rnd)]; end
    case 12
        if (rand()<.5), phi = [rnd; rnd; rnd]; else phi = [rnd; rnd; (1 - rnd)]; end
    case {13, 14}
        phi = rand(epsilon,1);
    case 16
        phi = [1; 1; 0];
    case 17
        phi = [1; 1; -1];
    case 18
        phi = [rnd; 1-rnd; rnd; 1-rnd];
end

%determine mask
switch scenario
    case 7
        mask = mask_1;
    case {8, 9, 10, 12, 16}
        if (rand()<.5), mask = mask_1; else mask = mask_2; end
    case 14
        mask = zeros(epsilon,1);
        mask(randi(epsilon))=1;
    case 18
        if (rand()<.5), mask = mask_3; else mask = mask_4; end
        
    otherwise
        mask = ones(epsilon,1); %default
end

if (scenario ~= 17)
    target = phi * (max_phi - min_phi) + min_phi; %adjust
else
    target = phi;
end

%estimate env cue
cue = (target - min_phi)*(max_cue - min_cue)/(max_phi - min_phi) + min_cue;
if (any(scenario == [0 1 15])), cue = 0 * target; end %no input
if (scenario == 17), cue = phi; end

end

