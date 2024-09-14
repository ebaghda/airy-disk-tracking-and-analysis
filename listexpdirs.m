function expDirectory = listexpdirs(dateFile)
% Takes files in format [mounted drive letter]:\[date]\[experiment] and
% returns a string array of the experiment folder names
dateDir = dir2(dateFile); % uses "dir2.m" function to ignore '.' and '..' folders
dateDir = dateDir(~cellfun('isempty',{dateDir.date})); %remove empty folders
isfolder = @(x) x == 1; %anonymous function to find folders
dateDir = dateDir(cellfun(isfolder,{dateDir.isdir})); %exclude non-folders
expDirs= [string({dateDir.folder})' repmat("\",size(string({dateDir.folder})')) string({dateDir.name})'];
expDirectory=erase(join(expDirs,2), " "); % list experiment paths
end