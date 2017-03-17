function [fullPath, slices] = getDCMslices(folderPath)
  currentDir = pwd;
  cd(folderPath)

  filecell = struct2cell([dir('*DCM*') dir('*MRDC*')]);

  fullPath = pwd;
  slices = filecell;
  cd(currentDir)
end
