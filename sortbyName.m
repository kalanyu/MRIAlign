function indices = sortByName(files)
  if isempty(files)
    indices = [];
  else
    extCell = cellfun(@(x) x(end - 1),cellfun(@(x) strsplit(x,'.'),files(1,:),'UniformOutput',false));
    [~, indices] = sort(cellfun(@(x) str2num(x), extCell),'ascend');
  end
end
