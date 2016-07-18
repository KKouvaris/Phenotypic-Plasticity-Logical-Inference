function [ folderDirectory ] = findFolder( scenario, N, t_max )

%Locate folder directory of the respective scenario

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


%Arguments validation
if ~exist('t_max','var') || isempty(t_max)
    t_max = 50000; %default
end
if ~exist('N','var') || isempty(N)
    warning('Error with population size');
    return;
end
if ~exist('scenario','var') || isempty(scenario)
    warning('Error with scenario number');
    return;
end

mainDirectory = 'D:\Data\Plasticity\Scenario_';

folderDirectory = [ mainDirectory num2str(scenario) '\' num2str(t_max) '\' num2str(N)];

end

