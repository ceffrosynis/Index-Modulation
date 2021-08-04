pkg load communications;
clc; clear all;

############# Parameters Section###########
CPL = 4;                        %Cyclic Prefix Length
NFFT = 64;                      %FFT size
modulation = 64;                %General modulation scheme
Eb = 1;                         %Power of bit
SNRdb = [0:5:30];               %SNR in db
############# Parameters Section###########


l = modulation / 4;             %Total size of subblock
M = modulation / l;             %Modulation scheme
symbolBits = log2(M);           %Bits per symbol
k = ones(1, l);
p = NFFT/l;                     %Number of subblocks

SNRdbSymbol = SNRdb + 10*log10(NFFT/NFFT) + 10*log10(64/80);          %OFDM symbol per noise ratio in db
SNRw = 10.^(SNRdbSymbol/10);                                          %OFDM symbol per noise ratio in Watt
Noise = 2*Eb./SNRw;                                                   %Power of Noise in Watt                                          

g1 = floor(log2(factorial(l)));       %Index bits  
g2 = l*log2(M);                                  %Data bits
g = g1 + g2;                                                        %Total bits for subblock                               
m = g * p;                                                          %Total bits for the FFT block
noblocks = 300;                                                      %Number of FFT blocks
N = m * noblocks;                                                   %Total number of input bits  

maxNumber = 2^g1;

input = randi([0, 1], 1, N);                                        %Input bits
if exist ('trellisDiagram.txt') == 2
  load('trellisDiagram.txt')
  'yes'
else  
  [stageIndexes, prevStageIndexes] = trellisGenerator (l);
endif  

symbols = generateSymbolGroups (modulation);

split = reshape(input, m, noblocks);

split = reshape(split, g, p, noblocks);

indexAdec1 = bin2dec(num2str(reshape(split(1:g1, :, :), g1, p*noblocks).'));

realIndexBits = reshape(split(1:g1, :, :), g1, p*noblocks).';

indexAdec = indexMapping(indexAdec1, l);

modeIndex = getIndexes (p*noblocks, l, indexAdec, k);

inputSymbols = reshape(split((g1+1):g, :, :), g2, noblocks*p)';

realData = inputSymbols;

inputSymbols = bin2dec(num2str(reshape(inputSymbols.', symbolBits, noblocks*p*g2/symbolBits).'));
inputSymbols = reshape(inputSymbols, p*noblocks, l);



inputSymbols = inputSymbols.';

modeGroup = cell(1,l);
ifftInput(1:NFFT*noblocks) = 0.';

for Mode = 1:l
  modeGroup(Mode) = inputSymbols(modeIndex{Mode}.');
  modeGroup(Mode) = symbols(Mode,modeGroup{Mode}+1);
  
  ifftInput(modeIndex{Mode}.') = modeGroup{Mode}.';
endfor

ifftInput = reshape(ifftInput, NFFT, noblocks).';

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

vars = {'split', 'ifftSignal', 'x', 'y', 'r', 'rt', 'symbolsA', 'symbolsB', 'indexA', 'indexB', 'inputSymbols', 'modeGroup', 'ifftInput'};
clear(vars{:})

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

  receivedSymbols = reshape(receivedSymbols, l, p*noblocks).';

  [finalIndex final] = ml(receivedSymbols, symbols, l, stageIndexes, prevStageIndexes);
  
  demodSymbols = [];
  
  for Mode = 1:l
    finalIndexA = finalIndex == Mode;
    demodSymbols(finalIndexA) = genqamdemod (final(finalIndexA), symbols(Mode,:));
  endfor

  finalModeIndex = indexDemapping(finalIndex, l);
  
  indexErrors = [indexErrors nnz(indexAdec1- finalModeIndex)];
  rightIndexData = find(~(indexAdec1- finalModeIndex));
  
  demodSymbols = reshape(demodSymbols, size(receivedSymbols,1), l);
  
  finalModeIndex = de2bi(finalModeIndex, g1, "left-msb");
  finalData = de2bi(demodSymbols, symbolBits , "left-msb");
  
  finalData = reshape(finalData.', l*symbolBits, size(receivedSymbols,1)).';
  finalModeIndex = reshape(finalModeIndex.', g1, size(receivedSymbols,1)).';
  
  indexBitsErrors= [indexBitsErrors nnz(finalModeIndex- realIndexBits)];
  
  rightIndexDataErrors = [rightIndexDataErrors nnz(finalData(rightIndexData, :)- realData(rightIndexData, :))];
  dataErrors = [dataErrors nnz(finalData- realData)];
  
  errors = [errors nnz([finalModeIndex finalData].'(:) - input.')];
  
endfor

dataMatSize = size(realData, 1) * size(realData, 2);
rightIndexDataErrors = rightIndexDataErrors / dataMatSize;

errors = errors/N;
%save('errors.txt', 'errors');

indexErrors = indexErrors/(size(indexAdec1, 1) * size(indexAdec1, 2));
%save('indexErrors.txt', 'indexErrors');

%save('rightIndexDataErrors.txt', 'rightIndexDataErrors');

semilogy(SNRdb,errors);
legend("16QAM 4-modes");








