function idxNumber = demap(indexArray, l)
  
  positions = 1:l;
  idxNumber = 0;
  noPermutations = factorial(l);
  val = noPermutations;
  
  for i = 1:l-1
    val = val / (l-i+1);
    idxNumber = idxNumber + (find(positions == indexArray(i)) - i) * val;
    positions(positions == indexArray(i)) = positions(i);
    positions(i) = indexArray(i);
  endfor  
end
