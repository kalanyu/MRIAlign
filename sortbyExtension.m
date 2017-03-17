function indices = sortByExtension(files)
  if isempty(files)
    indices = [];
  else
    extCell = cellfun(@(x) x(end),cellfun(@(x) strsplit(x,'.'),files(1,:),'UniformOutput',false));
    [~, indices] = sort(cellfun(@(x) str2num(x), extCell),'ascend');
  end
end
