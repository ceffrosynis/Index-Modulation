pkg load communications;
clc; clear all;
############# Parameters Section###########

CPL = 4;                %Cyclic Prefix Length
modulation = 16;
NFFT = 64;              %FFT size
l = 4;
k = 2;

symbols = qammod(0:(modulation-1), modulation);
indexSymbol = real(symbols) + 3 == 0 | real(symbols) + 1 == 0;


symbolBits = log2(modulation/2);        %Bits per symbol
MA = modulation / 2;                    %Modulation scheme A
MB = modulation / 2;                    %Modulation scheme B

############# Parameters Section###########

Eb = 1;                 %Power of bit
SNRdb = [0:10];           %SNR in db
SNRdbSymbol = SNRdb + 10*log10(NFFT/NFFT) + 10*log10(64/80);          %OFDM symbol per noise ratio in db
SNRw = 10.^(SNRdbSymbol/10);                                          %OFDM symbol per noise ratio in Watt
Noise = 2*Eb./SNRw;     %Power of Noise in Watt                                          

p = NFFT/l;                             %Number of subblocks

g1 = floor(log2(factorial(l)/(factorial(l-k)*factorial(k))));       %Index bits  

g2 = k*log2(MA) + (l-k)*log2(MB);                                   %Data bits

g = g1 + g2;                                                        %Total bits for subblock                               

m = g * p;                                                          %Total bits for the FFT block
noblocks = 50;                                                      %Number of FFT blocks
N = m * noblocks;                                                   %Total number of input bits  

maxNumber = 2^g1;

input = randi([0, 1], 1, N);                                        %Input bits

indeces = nchoosek([1:l], k);                                   
noindeces = nchoosek(l, k);

%%%%indexSymbol 

%Using the modulation scheme from the paper

symbols = [symbols(~indexSymbol); symbols(indexSymbol)];

split = reshape(input, m, noblocks);

split = reshape(split, g, p, noblocks);

indexAdec = bin2dec(num2str(reshape(split(1:g1, :, :), g1, p*noblocks).'));
indexAdec = indeces(indexAdec+1, :);

inputSymbols = reshape(split((g1+1):g, :, :), g2, noblocks*p)';

inputSymbols = bin2dec(num2str(reshape(inputSymbols.', symbolBits, noblocks*p*g2/symbolBits).'));
inputSymbols = reshape(inputSymbols, p*noblocks, l);

index = [1:l].*ones(p*noblocks, l);

c = indexAdec(:,1)-index == 0;
for i=2:k
  c = c | indexAdec(:,i)-index == 0;
endfor

indexA = c;
indexB = ~c;

inputSymbols = inputSymbols.';

%Splitting the MA symbols from MB
symbolsA = inputSymbols(indexA.');
symbolsA = symbols(2,symbolsA+1);
%symbolsA = reshape(symbolsA, k, p*noblocks).';

symbolsB = inputSymbols(indexB.');
symbolsB = symbols(1,symbolsB+1);
%symbolsB = reshape(symbolsB, l-k, p*noblocks).';

%Creating the FFT block
ifftInput(1:NFFT*noblocks) = 0.';
ifftInput(indexA.') = symbolsA.';
ifftInput(indexB.') = symbolsB.';

%ifftInput = reshape(ifftInput, l, p*noblocks).';

ifftInput = reshape(ifftInput, NFFT, noblocks).';

ifftSignal = ((NFFT/sqrt(NFFT)) * ifft(fftshift(ifftInput.'))).';

%Adding Cyclic Prefix
channelInput = ifftSignal(:, end-CPL+1:end);
channelInput = [channelInput ifftSignal];

totalLength = noblocks*(NFFT+CPL);

channelInput = reshape(channelInput.', 1, totalLength);


######### RICIAN ###############
Kdb = 10;  %Ratio of direct path and scattered paths
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

vars = {'input','split', 'ifftSignal', 'x', 'y', 'r', 'rt', 'ifftInput', 'symbolsA', 'symbolsB', 'indexA', 'indexB'};
clear(vars{:})

######### Clearing some space ##############


errors = [];

for noise = Noise
  demodSymbols = [];
  awgnnoise = sqrt(noise/2)*randn(1, totalLength) + j*sqrt(noise/2)*randn(1, totalLength);

  receivedSignal = h .* channelInput + awgnnoise;
  
  receivedSignal = receivedSignal ./ h;
  
  receivedSignal = reshape(receivedSignal.', NFFT+CPL, noblocks).';

%Removing cyclic prefix
  receivedSignal = receivedSignal(:, end-NFFT+1:end);

  receivedSymbols = (sqrt(NFFT)/NFFT)*fftshift(fft(receivedSignal.'));

  receivedSymbols = reshape(receivedSymbols, l, p*noblocks).';

  [c, i] = mlfun(receivedSymbols, l, k, p*noblocks, symbols);
  
  clear mlfun;
  
  cols = repmat([1:size(receivedSymbols,1)], k, 1);
  rows = indeces(i,:).';
  linearIndex = sub2ind(size(receivedSymbols.'), rows, cols);
  
  finalData = rows.';
  
  demodSymbolsA = receivedSymbols.'(linearIndex);
  demodSymbols(linearIndex) = genqamdemod (demodSymbolsA, symbols(2,:));
  
  cols = repmat([1:size(receivedSymbols,1)], l-k, 1);
  rows = repmat([1:l].', 1, size(rows, 2));
  rows(linearIndex) = [];
  rows = reshape(rows, l-k, size(cols, 2));
  linearIndex = sub2ind(size(receivedSymbols.'), rows, cols);
  
  demodSymbolsB = receivedSymbols.'(linearIndex);
  demodSymbols(linearIndex) = genqamdemod (demodSymbolsB, symbols(1,:));
  
  demodSymbols = reshape(demodSymbols.', l, size(receivedSymbols,1)).';
  
  finalData = [finalData demodSymbols];
  
  errors = [errors nnz(finalData-[indexAdec inputSymbols.'])];
  
endfor

theory = 0.5*erfc(sqrt(Eb./Noise));
semilogy(SNRdb,errors/(size(demodSymbols,1)*(size(demodSymbols,2)+k)), SNRdb, theory);
legend("8-PSK AWGN Simulation", "8-PSK AWGN Theory");








