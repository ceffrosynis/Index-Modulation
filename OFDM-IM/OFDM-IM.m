pkg load communications;
clc; clear all; close all;


####Arguments####
SNRdb = [0:10];     %SNR range
modulation = 16;    %Modulation scheme
N = 5000;          %Input length

NFFT = 60;     %Number of OFDM subcarriers  (64)
CPL = 16;      %Cyclic Prefix (16)

g = 10;  %Number of subblocks
m = 190;  %Size of info bits

k = 4;   %Number of activated subcarriers

Eb = 1;

######OFDM#######
DSC = k*g;     %Data subcarriers (64)
SNRdbSymbol = SNRdb + 10*log10(DSC/NFFT) + 10*log10(NFFT/(NFFT+CPL)); %Symbol per noise ratio
symbolbits = log2(modulation);

######INDEX-MODULATION#######
p = m/g;  %Number of bits in each subblock
n = NFFT/g; %OFDM-IM subblock length
p1 = k * symbolbits;  %M-ary signal
p2 = floor(log2(factorial(n)/(factorial(n-k)*factorial(k))));       %Index bits  

maxNumber = 2^p2;       %all possible combinations of carrier indexes
indeces = nchoosek([1:NFFT/g], k);
noindeces = nchoosek(NFFT/g, k);

N = N * m;

SNRw = 10.^(SNRdbSymbol/10);
Noise = 2*Eb./SNRw;

input = randi([0, 1], 1, N);

noblocks = N/m;         %Number of OFDM-IM blocks

%Splitting the input into a three dimensional array
%1st Dimension (Vertical) : Bits of each subblock
%2nd Dimension (Horizontal) : Splitting the subblocks
%3nd Dimension (z) : Splitting the OFDM-IM symbols

bitBlocks = reshape(input, m, noblocks);
bitBlocks = reshape(bitBlocks, p, g, noblocks);

activeSubcarriers = bin2dec(num2str(reshape(bitBlocks((p1+1):p, :, :), p2, noblocks*g)'));

realIndexBits = reshape(bitBlocks((p1+1):p, :, :), p2, noblocks*g);
realIndex = activeSubcarriers;


activeSubcarriers = indeces(activeSubcarriers+1, :);    %Mapping the p2 bits into active subcarrier indeces

realData = reshape(bitBlocks(1:p1, :, :), p1, noblocks*g).';   

inputSymbols = bin2dec(num2str(reshape(bitBlocks(1:p1, :, :), symbolbits, noblocks*g * k)'));

%Encoding M-ary modulated data
complexInputSymbols = qammod(inputSymbols, modulation);

offset = [0:n:noblocks*NFFT-n]';
activeSubcarriers = activeSubcarriers + offset;
%Initializing the fft blocks
ifftInput = zeros(noblocks, NFFT);

%Placing the data in the appropriate positions
ifftInput(activeSubcarriers) = complexInputSymbols;

%Applaying the ifft
ifftSignal = ((NFFT/sqrt(DSC)) * ifft(fftshift(ifftInput.'))).';
%Adding Cyclic Prefix
SignalWithCyclicPrefix = ifftSignal(:, end-CPL+1:end);
SignalWithCyclicPrefix = [SignalWithCyclicPrefix ifftSignal];

%Serializing the data
totalLength = noblocks*(NFFT+CPL);
channelInput = reshape(SignalWithCyclicPrefix.', 1, totalLength);

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

errors = []
rightIndexDataErrors = [];

for noise = Noise
  
  %AWGN and Rician channel fading
  awgnnoise = sqrt(noise/2)*randn(1, totalLength) + j*sqrt(noise/2)*randn(1, totalLength);
  data = h .* channelInput + awgnnoise;
  
  %Removing the channel fading (We assume that the receiver knows the channel coefficients and the CHANNEL TAPS < CP)
  data = data ./ h;

  receivedSignal = reshape(data.', NFFT+CPL, noblocks).';
  
  %Removing Cyclic Prefix
  receivedSignal = receivedSignal(:, end-NFFT+1:end);
  
  %Applying fft
  receivedSymbols = (sqrt(DSC)/NFFT)*fftshift(fft(receivedSignal.'));
  
  receivedSymbols = reshape(receivedSymbols.', n, noblocks*g);
  
  %Sorting the data by power in order to find the most probable active carriers
  [sorted, index] = sort(abs(receivedSymbols), "descend");
  blockIndex = index(1:k, :);
  [sorted, index] = sort(blockIndex, "ascend");
  
  pos = zeros(maxNumber, noblocks*g);
  %Finding the indeces of active subcarriers
  for i = 1:k
    
    a = sorted(i,:) - indeces (1:maxNumber,i);
    pos = pos | a;
    
  endfor
  
  [pos ic] = min(pos);
  
  d1 = dec2bin(ic-1);
  d1 = str2num(d1(:));
  d1 = reshape(d1, noblocks*g, p2);
  
  newIndex = sorted.' + offset;
  %Retrieving the M-ary data
  receivedSymbols = receivedSymbols(newIndex);
  %Demodulate M-ary data
  data = qamdemod(receivedSymbols, modulation);
  
  data = reshape(data, g*k, noblocks);
  
  finalData = de2bi(data, symbolbits, "left-msb");
  finalData= reshape(finalData.', p1, noblocks* g).';
  
  %Concatenating the active indeces with the data
  
  rightIndexData = find(~([ic - 1].' - realIndex));
  
  data = [finalData d1];
  data = reshape(data.', 1, N);

  errors = [errors nnz(input-data)];

endfor

dataMatSize = size(realData, 1) * size(realData, 2);

errors = errors/N;

semilogy(SNRdb,errors);

legend("BPSK OFDM-IM Simulation");

