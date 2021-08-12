pkg load communications;
clc; clear all;

addpath("./../lib/IndexMapper");           % Library path

############# Parameters Section###########

CPL = 4;                              %Cyclic Prefix Length
modulation = 16;                      %General modulation order
l = modulation / 4;                   %Total size of sub-block
M = modulation / l;                   %Modulation order corresponding to each constellation signal
u = 2;                                %Number of sub-sub-blocks
p = 2;                                %Number of sub-blocks per FFT block
Eb = 1;                               %Power of bit in Watt
SNRdb = [1:20];                         %SNR range of interest in db
noblocks = 100;                         %Number of blocks per SI-MM-OFDM-IM symbol
############# Parameters Section###########

NFFT = p * l;                         %FFT size
DSC = NFFT; 			
symbolBits = log2(M);                 %Bits per symbol
k = ones(1, l);

SNRdbSymbol = SNRdb + 10*log10(DSC/NFFT) + 10*log10(NFFT/(NFFT+CPL));          %OFDM symbol per noise ratio in db
SNRw = 10.^(SNRdbSymbol/10);                                                   %OFDM symbol per noise ratio in Watt
Noise = 2*Eb./SNRw;                                                            %Power of Noise in Watt                                          

g1 = floor(log2(factorial(u)));                                     %Sub-block Index Bits
g2 = u * floor(log2(factorial(l/u)));                               %Mode Index Bits
g3 = l*log2(M);                                                     %Data Bits
g = g1 + g2 + g3;                                                   %Total bits for subblock                               

m = g * p;                                                          %Total bits of the FFT block

N = m * noblocks;                                                   %Total number of input bits  

input = randi([0, 1], 1, N);                                        %Input bits

if and(exist('trellisDiagram.txt'), l/u==64)
  load('trellisDiagram.txt')
  subSubBlockIndexes = stageIndexes;
  subSubBlockprevIndexes = prevStageIndexes;
  [subBlockIndexes, subBlockprevIndexes] = trellisGenerator (u);
elseif and(exist('trellisDiagram.txt'), u==64)
  load('trellisDiagram.txt')
  subBlockIndexes = stageIndexes;
  subBlockprevIndexes = prevStageIndexes;
  [subSubBlockIndexes, subSubBlockprevIndexes] = trellisGenerator (u);
else  
  [subSubBlockIndexes, subSubBlockprevIndexes] = trellisGenerator (l/u);
  [subBlockIndexes, subBlockprevIndexes] = trellisGenerator (u);
endif

symbols = generateSymbolGroups (modulation);

split = reshape(input, m, noblocks);

split = reshape(split, g, p, noblocks);

blockIndexBits = bin2dec(num2str(reshape(split(1:g1, :, :), g1, p*noblocks).'));

blockIndexDec = bin2dec(num2str(reshape(split(1:g1, :, :), g1, p*noblocks).'));

realIndexBits = reshape(split(1:g1, :, :), g1, p*noblocks).';

blockIndexDec = indexMapping(blockIndexDec, u);
if length(blockIndexDec) == 0
  blockIndexDec= 1;
endif

indexBits = reshape(split((g1+1):(g1+g2), :, :), g2, noblocks*p)';

indexDec = bin2dec(num2str(reshape(indexBits.', g2/u, noblocks*p*u).'));

indexDec = indexMapping(indexDec, l/u);

indeces = (repmat(blockIndexDec.'(:), 1, l/u) - 1) * (l/u) + indexDec;

dataBits = reshape(split((g1+g2+1):(g), :, :), g3, noblocks*p)';

dataDec = bin2dec(num2str(reshape(dataBits.', g3/l, noblocks*p*l).'));

dataDec = reshape(dataDec.', l/u, u*p*noblocks).';

indexSymbols = sub2ind (size(symbols), indeces, dataDec+1);

ifftInput = symbols(indexSymbols);

ifftInput = reshape(ifftInput.', NFFT, noblocks).';

ifftSignal = ((NFFT/sqrt(NFFT)) * ifft(fftshift(ifftInput.'))).';

%Adding Cyclic Prefix
channelInput = ifftSignal(:, end-CPL+1:end);
channelInput = [channelInput ifftSignal];

totalLength = noblocks*(NFFT+CPL);

channelInput = reshape(channelInput.', 1, totalLength);


######### RICIAN ###############
Kdb = 2;  %Ratio of direct path and scattered paths
Mi = 1;   %Interpolation factor

K = 10 ^ (Kdb/10);
const = 1/(2*(K+1));
x = randn(1,totalLength);
y = randn(1,totalLength);
r = sqrt(const*((x + sqrt(2*K)).^2 + y.^2));
rt = zeros(1,Mi*length(r));
ki = 1;

for i=1:length(r)
  rt(ki:i*Mi) = r(i);
  ki = ki+Mi;
end
h = rt;

##################################

######### Clearing some space ##############



######### Clearing some space ##############


errors = [];
indexErrors = [];
indexBitsErrors = [];
dataErrors = [];
rightIndexDataErrors = [];



for noise = Noise
  noise
  demodSymbols = [];
  awgnnoise = sqrt(noise/2)*randn(1, totalLength) + j*sqrt(noise/2)*randn(1, totalLength);

  receivedSignal = h .* channelInput + awgnnoise;
  
  receivedSignal = receivedSignal ./ h;
  
  receivedSignal = reshape(receivedSignal.', NFFT+CPL, noblocks).';
  
  receivedSignal = receivedSignal(:, end-NFFT+1:end);

  receivedSymbols = (sqrt(NFFT)/NFFT)*fftshift(fft(receivedSignal.'));
  
  clear receivedSignal;

  receivedSymbols = reshape(receivedSymbols .',l/u, p*noblocks*u ).';
  
  baseIndex = 1:l/u;
  
  matrixSymbols = [];
  matrixModes = [];
  
  for mode = 1:u
    mode
    matrixIndex = (mode - 1)*l/u + baseIndex;
    [finalIndex final] = ml(receivedSymbols, symbols(matrixIndex,:), l/u, subSubBlockIndexes, subSubBlockprevIndexes);
    matrixSymbols = cat(3, matrixSymbols, final);
    matrixModes = cat(3, matrixModes, finalIndex);
  endfor  
  'mlmode'
  [finalIndex ] = mlMode (receivedSymbols, matrixSymbols, u, l, subBlockIndexes, subBlockprevIndexes);
  
  matrixIndeces = sub2ind (size(matrixModes), repmat([1:size(receivedSymbols, 1)].', 1, l/u), repmat([1:l/u], size(receivedSymbols, 1), 1), repmat(finalIndex.'(:), 1, l/u));
  finalData = matrixSymbols(matrixIndeces);
  
  finalData = genqamdemod (finalData, symbols.'(:));
  finalData= mod(finalData, M);
  
  finalDataBits = de2bi (finalData.', g3/l, "left-msb");
  finalDataBits = reshape(finalDataBits.', g3, p*noblocks).';
  
  finalModes = matrixModes(matrixIndeces);
  
  finalModeIndexBits = indexDemapping(finalModes, l/u);
  
  finalModeIndexBits = de2bi(finalModeIndexBits, g2/u, "left-msb");
  finalModeIndexBits = reshape(finalModeIndexBits.', g2, p*noblocks).';
  
  finalBlockModeIndex = indexDemapping(finalIndex, u);
  
  finalBlockModebits = de2bi(finalBlockModeIndex, g1, "left-msb");
  
  data = [finalBlockModebits finalModeIndexBits finalDataBits];
  
  data = data.'(:).';
  
  errors = [errors nnz(input- data)]
  
endfor

errors = errors/N;

semilogy(SNRdb,errors);
legend("16QAM 4-modes");



