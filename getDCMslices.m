function [fullPath, slices] = getDCMslices(folderPath)
  currentDir = pwd;
  cd(folderPath)

  filecell = struct2cell([dir('*DCM*') dir('*MRDC*')]);

  fullPath = pwd;
  nonhiddenIndices = find(arrayfun(@(x) ~strcmp(x{:}(1), '.'), filecell(1,:)));
  slices = filecell(:,nonhiddenIndices);
  
  cd(currentDir)
end
