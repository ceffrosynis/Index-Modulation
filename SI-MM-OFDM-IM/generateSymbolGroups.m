% This function generates the constellation symbol groups based on the analysis of the Constellation Design.

function groupedSymbols = generateSymbolGroups (modulation)
  
  allPossibleQAMSchemes = [4 16 64 256 1024];
  
  diffValues = allPossibleQAMSchemes - modulation;
  
  [maxVal idx] = max(diffValues(diffValues <= 0));
  M = allPossibleQAMSchemes(idx);
  
  symbols = qammod(0:allPossibleQAMSchemes(idx+1)-1, allPossibleQAMSchemes(idx+1));
  
  groupedSymbols = [];
  
  realPart = 1:2:2*sqrt(M)/2;
  realPart = [-realPart(end:-1:1).' realPart.'];
  
  imagPart = 1:2:modulation/sqrt(M);
  imagPart = -1.*[-imagPart(end:-1:1).' imagPart.'];
  
  realPoints = size(realPart, 1);
  imagPoints = size(imagPart, 1);
  
  for iPart = 1:imagPoints
    for rPart = 1:realPoints
      tmp = (real(symbols) == realPart(rPart, 1) | real(symbols) == realPart(rPart, 2)) & (imag(symbols) == imagPart(iPart, 1) | imag(symbols) == imagPart(iPart, 2));
      groupedSymbols = [groupedSymbols; symbols(tmp)];
    endfor
  endfor
  
  
endfunction
