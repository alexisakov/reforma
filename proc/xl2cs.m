function [ output_args ] = xl2cv( file,sheet,exportname )
%XL2CV Summary of this function goes here
%   Detailed explanation goes here

[pathstr,~,~] = fileparts(file);
[~, ~, rawData] = xlsread(file,sheet);
rawData{1,1}='';
 fid = fopen(fullfile(pathstr, exportname), 'w') ;
 for ix = 1:size(rawData,1)
     fprintf(fid, '%s,', rawData{ix,1:end-1}) ;
     fprintf(fid, '%s\n', rawData{ix,end}) ;
 end
 fclose(fid) ;
end

