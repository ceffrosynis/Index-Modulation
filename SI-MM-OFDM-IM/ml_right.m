%Performs the Maximum Likelihood Detection over the received SI-MM-OFDM-IM symbols (based on the Viterbi - like algorithm).

function [finalIndex final] = ml (a, b, l, stageIndexes, prevStageIndexes)

  noblocks = size(a, 1);

  row = size(a, 1);
  column = size(a, 2);

  symbols = [];
  symbolArray = [];
  stageMin = cell(l, 1);
  stageMinIndex = cell(l, 1);
  finalState = cell(l, 1);
  linIndex = cell(l, 1);
  linPrevStage = cell(l, 1);

  for i = 1:size(b, 1)
    q = [];
    for u = 1:size(b, 2)
      q = cat(3, q, abs(a-b(i,u)));
    endfor
    [val, symbolIndex] = min(q, [], 3);
    symbols = cat(2, symbols, reshape(min(q, [], 3).', [column 1 row]));
    symbolArray = cat(2, symbolArray, reshape(b(i,symbolIndex.'), [column 1 row]));
  endfor
  
  stageMin(1)=permute(symbols(1, stageIndexes{1}, :), [2 1 3]);


  for stage = 2:l-1
    symbols2 = permute(symbols(stage, stageIndexes{stage}, :), [2 1 3]) + stageMin{stage-1}(prevStageIndexes{stage},:,:);
    symbols2 = reshape(symbols2, [size(stageIndexes{stage}) noblocks]);
    [stageMin(stage) stageMinIndex(stage)] = min(symbols2, [], 2);
  endfor

  symbols4=permute(symbols(l, stageIndexes{l}, :), [2 1 3])+stageMin{l-1}(prevStageIndexes{l},:,:);
  [stageMin(l) stageMinIndex(l)] = min(symbols4, [], 1);


  linIndex4 = sub2ind(size(stageMinIndex{l}), ones(noblocks, 1), ones(noblocks, 1), [1:noblocks].');
  finalState(l) = permute(stageMinIndex{l}(linIndex4), [3 1 2]);
  stageMinIndex(l) = finalState{l};

  for stage = l-1:-1:2
    linIndex(stage) = sub2ind(size(stageMinIndex{stage}), finalState{stage+1}, ones(noblocks, 1), [1:noblocks].');
    stageMinIndex(stage) = stageMinIndex{stage}(linIndex{stage});

    linPrevStage(stage) = sub2ind(size(prevStageIndexes{stage}), finalState{stage+1}, stageMinIndex{stage});
    finalState(stage) = prevStageIndexes{stage}(linPrevStage{stage});
  endfor

  for stage = 2:l-1    
    stageMinIndex(stage) = sub2ind(size(stageIndexes{stage}), finalState{stage+1}, stageMinIndex{stage});
  endfor
  
  finalIndex = [];
  stageMinIndex(1) = finalState{2};

  for stage = 1:l    
    finalIndex = [finalIndex stageIndexes{stage}(stageMinIndex{stage})(:)];
  endfor
  
  linearFinal = sub2ind(size(symbolArray), repmat(1:l, noblocks, 1), finalIndex, repmat([1:noblocks].', 1, l));

  final = [symbolArray(linearFinal)];

endfunction