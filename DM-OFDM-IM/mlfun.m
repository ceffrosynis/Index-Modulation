function [minVal, idx] = mlfun (a, l, k, nosubs, b1)
  row1 = reshape(repmat([1:nosubs].', 1, k).', nosubs*k, 1);
  row2 = reshape(repmat([1:nosubs].', 1, l-k).', nosubs*(l-k), 1);
  
  vars = {'c1','c2'};

  columns = nchoosek([1:l], k);
  nocolumns = nchoosek(l, k);

  ######## Finding all possible combinations of symbolsA #######
  C = cell(1,k);
  [C{:}] = ndgrid(b1(2,:));
  C = cellfun(@(a)a(:),C,'Uni',0);
  MA = [C{:}];
  ######## Finding all possible combinations of symbolsA #######
  
  ######## Finding all possible combinations of symbolsB #######
  C = cell(1,l-k);
  [C{:}] = ndgrid(b1(1,:));
  C = cellfun(@(a)a(:),C,'Uni',0);
  MB = [C{:}];
  ######## Finding all possible combinations of symbolsB #######
  
  dataSize = size(a, 1);     %size of data
  modSizeA = size(MA, 1);    %size of modulation scheme
  modSizeB = size(MB, 1);    %size of modulation scheme
  
  minVal = [];
  for row = 1:nocolumns
    ######## Fin #######
    rowsA = repmat(columns(row, :), 1, nosubs).';
    linearIndexA = sub2ind(size(a), row1, rowsA);
    linearIndexA = reshape(linearIndexA, k, nosubs).';
    otherColumns = 1:l;
    otherColumns(columns(row,:)) = [];
    rowsB = repmat(otherColumns, 1, nosubs).';
    linearIndexB = sub2ind(size(a), row2, rowsB);
    linearIndexB = reshape(linearIndexB, l-k, nosubs).';
    ######## Finding all possible combinations of symbolsA #######
    %% Assuming that the l equals k
    c1 = [];
    c2 = [];
    for w = 1:l-k
        c1 = [c1 abs(a(linearIndexB)(:,w) - MB(:,w).')(:).^2];
    endfor
    
    for w = 1:k          
        c2 = [c2 abs(a(linearIndexA)(:,w) - MA(:,w).')(:).^2];
    endfor
    
    c1 = sum(c1, 2);
    c2 = sum(c2, 2);
    
    %tic
    c1 = cellslices(c1, 1:dataSize:size(c1, 1), dataSize:dataSize:size(c1, 1));
    c2 = cellslices(c2, 1:dataSize:size(c2, 1), dataSize:dataSize:size(c2, 1));
    
    c1 = repmat(c1, modSizeA, 1);
    c2 = repmat(c2, modSizeB, 1).';
    %for cellElem = 1:modSize*modSize
      %C = sum([c1{cellElem} c2{cellElem}], 2);
    %endfor  
    
    %C = arrayfun(@(k) sum([c1{k} c2{k}], 2), 1:modSize*modSize, 'un', 0);
    %c1 = reshape(c1, modSizeA, modSizeA);
    %c2 = reshape(c2, modSizeA, modSizeA);
    
    %C = reshape(C, modSize, modSize);
    c1 = cellfun(@plus, c1, c2, 'Uni', 0);
    c1 = horzcat(c1{});
    
    %toc
    minVal = [minVal min(c1, [], 2)];
    clear(vars{:})
    
  endfor
  
  [minVal, idx] = min(minVal, [], 2);
  
endfunction
