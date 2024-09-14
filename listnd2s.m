function [nd2Files, nd2Count, fileFullPaths] = listnd2s(expDirs,eNum)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%expDir =dir2(expDirs(eNum));
files = dir2(expDirs(eNum));
fileDirs= [string({files.folder})' repmat("\",size(string({files.folder})')) string({files.name})'];
fileFullPaths=erase(join(fileDirs), " ");
files=struct2table(files);
nd2Count = find(endsWith(files.name, ".nd2"));
nd2Files=files(endsWith(files.name, ".nd2"),:);

if numel(nd2Count)>1
    disp("Multiple .nd2 files found in the current experiment")
    head(nd2Files)
end



end