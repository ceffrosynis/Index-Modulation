function modeIndex = getIndexes (dataSize, blockLength, indexArray, k)
  
  
  cumulativeSum = cumsum([0 k], 2);
  index = [1:blockLength].*ones(dataSize, blockLength);
  
  modeIndex = cell(1, blockLength);
  
  for Mode = 1:blockLength
    modeIndex(Mode) = indexArray-index(:,cumulativeSum(Mode)+1) == 0;
    for i = cumulativeSum(Mode)+2:cumulativeSum(Mode+1)
      modeIndex(Mode) = modeIndex{Mode} | indexArray(:,i)-index == 0;
    endfor
  endfor
  
endfunction
