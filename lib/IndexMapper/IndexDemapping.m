%Takes as an input each mode's index pattern and returns the index number.

function decimalIndex = indexDemapping(indexArray, l)
  
  curLength = l:-1:1;
  
  noPermutations = factorial(l);
  decimalIndex = zeros(size(indexArray), 1);
  val = noPermutations;
  block = zeros(size(indexArray));
  positions = repmat([1:l].', 1, size(indexArray, 1));
  
  for i = 1:l-1
    val = val / curLength(i);
    
    linearIndex1 = sub2ind(size(positions), indexArray(:,i).', [1:size(positions, 2)]);
    
    block(:,i) = positions(linearIndex1);

    positions(positions == i) = positions(linearIndex1);
    
    positions(linearIndex1) = i;
    
    decimalIndex = decimalIndex + val .* (block(:,i) - i);
  endfor
endfunction
