function [labels, pos] = readOptoPositions( filePath )
%READOPTOPOSITIONS Summary of this function goes here
%   Detailed explanation goes here
  filePath
  fid = fopen(filePath);
  contents = textscan(fid, '%f %f %f %f %f %f %f %f', 'Delimiter', ',');
  % celldisp(contents);
  output_args = cell2mat(contents);
  labels = output_args(:,1);
  pos = output_args(:,2:4);
end
