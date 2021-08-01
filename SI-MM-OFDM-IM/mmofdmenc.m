pkg load communications;
clc; clear all;

addpath("./../lib/IndexMapper");           % Library path

############# Parameters Section###########

CPL = 4;                          %Cyclic Prefix Length
modulation = 16;                  %General modulation scheme
l = modulation / 4;               %Total size of subblock
M = modulation / l;               %Modulation scheme
u = 2;                            %Number of sub- sub-blocks per group
p = 1;                            %Number of groups per FFT block
Eb = 1;                           %Power of bit in Watt
SNRdb = [0:5:30];                 %SNR range of interest in db
noblocks = 100;                     %Number of FFT blocks per SI-MM-OFDM-IM symbol
############# Parameters Section###########

NFFT = p * l;                     %FFT block size
DSC = NFFT;                       %Active OFDM subcarriers OFDM per sub-block

symbolBits = log2(M);        %Bits per symbol

k = ones(1, l);       

SNRdbSymbol = SNRdb + 10*log10(NFFT/NFFT) + 10*log10(64/80);          %OFDM symbol per noise ratio in db
SNRw = 10.^(SNRdbSymbol/10);                                          %OFDM symbol per noise ratio in Watt
Noise = 2*Eb./SNRw;     %Power of Noise in Watt                                          

g1 = floor(log2(factorial(u)));                                     %Sub-block Index Bits
g2 = u * floor(log2(factorial(l/u)));                                  %Mode Index Bits
g3 = l*log2(M);                                                  %Data Bits
g = g1 + g2 + g3;                                                        %Total bits for subblock                               

m = g * p;                                                          %Total bits for the FFT block
                                                     %Number of FFT blocks
N = m * noblocks;                                                   %Total number of input bits  

input = randi([0, 1], 1, N);                                        %Input bits

%% In case we have a trellis diagram with noStates == 16, we prefer to load a precomputed
%% diagram due to the excessive computational complexity.
if exist ('trellisDiagram.txt')  
  load('trellisDiagram.txt')
  if l/u == 16
    stageIndexesModes=stageIndexes;
    prevStageIndexesModes = prevStageIndexes;
    [stageIndexesGroups, prevStageIndexesGroups] = trellisGenerator (u);
  elseif u == 16
    stageIndexesGroups=stageIndexes;
    prevStageIndexesGroups = prevStageIndexes;
    [stageIndexesModes, prevStageIndexesModes] = trellisGenerator (l/u);
  else
    [stageIndexesModes, prevStageIndexesModes] = trellisGenerator (l/u);
    [stageIndexesGroups, prevStageIndexesGroups] = trellisGenerator (u);
  endif
else  
  [stageIndexesModes, prevStageIndexesModes] = trellisGenerator (l/u);
  [stageIndexesGroups, prevStageIndexesGroups] = trellisGenerator (u);
endif  


symbols = generateSymbolGroups (modulation);

inputBits = reshape(input, m, noblocks);
inputBits = reshape(inputBits, g, p, noblocks);

blockIndexDec = bin2dec(num2str(reshape(inputBits(1:g1, :, :), g1, p*noblocks).'));

realIndexBits = blockIndexDec;

blockIndexDec = IndexMapping(blockIndexDec, u);
if length(blockIndexDec) == 0
  blockIndexDec= 1;
endif

indexBits = reshape(inputBits((g1+1):(g1+g2), :, :), g2, noblocks*p)';

indexBits = bin2dec(num2str(reshape(indexBits.', g2/u, noblocks*p*u).'));

indexDec = IndexMapping(indexBits, l/u);

indeces = (repmat(blockIndexDec.'(:), 1, l/u) - 1) * (l/u) + indexDec;

dataBits = reshape(inputBits((g1+g2+1):(g), :, :), g3, noblocks*p)';

dataDec = bin2dec(num2str(reshape(dataBits.', g3/l, noblocks*p*l).'));

dataDec = reshape(dataDec.', l/u, u*p*noblocks).';

encryptedInputSymbols = dataDec + indexBits+ repmat(realIndexBits, u, 1);
encryptedInputSymbols = mod(encryptedInputSymbols, M);

indexSymbols = sub2ind (size(symbols), indeces, encryptedInputSymbols+1);

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
    [finalIndex final] = ml(receivedSymbols, symbols(matrixIndex,:), l/u, stageIndexesModes, prevStageIndexesModes);
    matrixSymbols = cat(3, matrixSymbols, final);
    matrixModes = cat(3, matrixModes, finalIndex);
  endfor  
  'mlmode'
  [finalIndex ] = mlMode (receivedSymbols, matrixSymbols, u, l, stageIndexesGroups, prevStageIndexesGroups);
  
  matrixIndeces = sub2ind (size(matrixModes), repmat([1:size(receivedSymbols, 1)].', 1, l/u), repmat([1:l/u], size(receivedSymbols, 1), 1), repmat(finalIndex.'(:), 1, l/u));
  finalData = matrixSymbols(matrixIndeces);
  
  finalData = genqamdemod (finalData, symbols.'(:));
  finalData = reshape(finalData.', l/u, u*p*noblocks).';
  
  finalModes = matrixModes(matrixIndeces);
  
  finalModeIndexBits = IndexDemapping(finalModes, l/u);
  
  finalModeKey = finalModeIndexBits;
  
  finalModeIndexBits = de2bi(finalModeIndexBits, g2/u, "left-msb");
  finalModeIndexBits = reshape(finalModeIndexBits.', g2, p*noblocks).';
  
  finalBlockModeIndex = IndexDemapping(finalIndex, u);
  
  finalBlockModeKey = finalBlockModeIndex;
  
  unencryptedSymbols = finalData - finalModeKey - repmat(finalBlockModeKey, u, 1);
  unencryptedSymbols = mod(unencryptedSymbols, M);
  
  finalData= mod(unencryptedSymbols, M);
  
  finalDataBits = de2bi (finalData.', g3/l, "left-msb");
  finalDataBits = reshape(finalDataBits.', g3, p*noblocks).';
  
  finalBlockModebits = de2bi(finalBlockModeIndex, g1, "left-msb");
  
  data = [finalBlockModebits finalModeIndexBits finalDataBits];
  
  data = data.'(:).';
  
  errors = [errors nnz(input- data)]
  
endfor

errors = errors/N;

semilogy(SNRdb,errors);
legend("16QAM 4-modes");



