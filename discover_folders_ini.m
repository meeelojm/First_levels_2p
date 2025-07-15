function [paths, fpss] = discover_folders_ini(root, defaultFPS)
%DISCOVER_FOLDERS_INI Locate subfolders containing both a .ini file and results.mat
%   [paths, fpss] = discover_folders_ini(root, defaultFPS)
    dirList = strsplit(genpath(root), pathsep);
    paths   = {};
    fpss    = [];
    for i = 1:numel(dirList)
        thisDir = dirList{i};
        if isempty(thisDir), continue; end
        % Look for at least one .ini file
        iniFiles = dir(fullfile(thisDir, '*.ini'));
        % Check for results.mat
        hasResults = isfile(fullfile(thisDir, 'results.mat'));
        if ~isempty(iniFiles) && hasResults
            paths{end+1} = thisDir;                %#ok<AGROW>
            fpss(end+1)  = read_fps(thisDir, defaultFPS);  %#ok<AGROW>
        end
    end
end
