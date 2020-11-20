%Generates the trellis diagram.

function [stageIndexes, prevStageIndexes] = trellisGenerator (l)
  
  stageIndexes = cell (l,1);
  prevStageIndexes = cell (l,1);
  
  stageIndexes(1) = flip(nchoosek([1:l], 1), 2);
  
  for stage = 2:l
    prevStates = [];
    stageIndexes(stage) = flip(nchoosek([1:l], stage), 2);
    for step = 1:size(stageIndexes{stage}, 1)
      prevStates = [prevStates; find(min(ismember(stageIndexes{stage-1}, stageIndexes{stage}(step, :)), [], 2)).'];
    endfor
    prevStageIndexes(stage) = prevStates;
  endfor
  
  
endfunction
