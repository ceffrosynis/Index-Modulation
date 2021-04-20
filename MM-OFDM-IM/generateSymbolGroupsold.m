function groupedSymbols = generateSymbolGroups (modulation)
  
  symbols = qammod(0:(modulation-1), modulation);
  
  groupedSymbols = [];
#{  
  realPart = 1:2:sqrt(modulation);
  realPart = [-realPart(end:-1:1).' realPart.'];
  
  imagPart = -realPart;
  #}
  
  realPart = 1:2:7;
  realPart = [-realPart(end:-1:1).' realPart.'];
  
  imagPart = 1:2:3;
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
