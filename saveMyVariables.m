function [] = saveMyVariables( scenario, saveObj ) %#ok<*INUSD>

%Arguments validation
if nargin < 2, return; end

Directory = findFolder(scenario, saveObj.N, saveObj.t_max);

if ~exist(Directory, 'dir')
  mkdir(Directory);
end

preFix = '\replicate_';

%find next entry
dirList = dir([Directory '\*.mat']);
numList = cell(length(dirList));

for id = 1:length(dirList)
    % Get the file name (minus the extension)
    [~, f] = fileparts(dirList(id).name);
    if ~isempty(strfind(f,'Pop'))
        findToken = strfind(f,'_');
        numList{id} = str2double(f(findToken(1)+1:findToken(2)-1));
    end
end

if isempty(dirList), entryNumber = 1; else entryNumber = max([numList{:}]) + 1; end

%save
if (isfield(saveObj,'Pop')), save([Directory preFix num2str(entryNumber) '_Pop.mat'],'-struct','saveObj','Pop'); end

end

