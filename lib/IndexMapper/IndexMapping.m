%Takes as an input the index number (index bits) and returns the index pattern.

function block =  indexMapping (index, l)
  
  length = l:-1:1
  noPermutations = factorial(l);
  valueRange = noPermutations;
  block = repmat(1:l, size(index, 1), 1);
  blockSize = size(block);
  
  for i = 1:l-1
    valueRange = valueRange / length(i);
    curIndex = floor(index ./ valueRange);
    linearIndex = sub2ind(blockSize, [1:blockSize(1)].', i+curIndex);
    curElement = block(:,i);
    block(:,i) = block(linearIndex);
    block(linearIndex) = curElement;
    index = mod(index, valueRange);
  endfor
  
endfunction
