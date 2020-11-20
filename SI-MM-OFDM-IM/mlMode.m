%Perfroms the Maximul Likelihood Detection over all possible constellations sets.

function [finalIndex ] = mlMode (a, b, u, l, stageIndexes, prevStageIndexes)

  noblocks = size(a, 1)/u;

  matrix = abs(a - b);
  blockSum = sum(matrix, 2);
  
  symbols = reshape(blockSum, u, noblocks, u);
  symbols = permute(symbols, [1 3 2]);

  row = size(a, 1);
  column = size(a, 2);

  l=u;
  
  stageMin = cell(l, 1);
  stageMinIndex = cell(l, 1);
  finalState = cell(l, 1);
  linIndex = cell(l, 1);
  linPrevStage = cell(l, 1);
  
  
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

endfunction