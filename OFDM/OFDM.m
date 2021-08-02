###### Initialization  ########
pkg load communications;
clc; clear all; close all;


#### Arguments ####
SNRdb = [0:5:30];	    ####### SNR Range (to db) ########
modulation = 4;		    ####### Modulation Scheme ########
noFFTBlocks = 1000;		####### Number of transmitting FFT Blocks #######
Rb = 10000;		        ####### Transmitting Bit Rate #######
oversampling = 1;	    ####### Oversampling at Transmitter/Receiver ########

###### OFDM #######
NFFT = 64;     		%FFT Size
CPL = 16;      		%Cyclic Prefix
DSC = 64;      		%Data (active) subcarriers
SNRdbSymbol = SNRdb + 10*log10(DSC/NFFT) + 10*log10(NFFT/(NFFT+CPL)); 	%Symbol per noise ratio

Eb = 1;			%The energy of one bit
symbolbits = log2(modulation);

N = NoFFTBlocks * symbolbits * NFFT;

SNRw = 10.^(SNRdbSymbol/10);
Noise = 2 * Eb ./ SNRw;

######## Normalized oversampling filters ########
txfilter = sqrt(Eb/oversampling)*ones(1, oversampling);
rxfilter = sqrt(Eb/oversampling)*ones(1, oversampling);

input = randi([0, 1], 1, N);

###### Data Modulation #######
inputSymbols = reshape(input, symbolbits, N/symbolbits).';
inputSymbols = bin2dec(num2str(inputSymbols));
complexInputSymbols = qammod(inputSymbols, modulation);  

nobins = length(complexInputSymbols)/NFFT;		      ## Number of FFT blocks
ifftInput = reshape(complexInputSymbols, NFFT, nobins).';

###### FFT Operation + CP Insertion ######
ifftSignal = ((NFFT/sqrt(DSC)) * ifft(fftshift(ifftInput.'))).';
SignalWithCyclicPrefix = ifftSignal(:, end-CPL+1:end);
SignalWithCyclicPrefix = [SignalWithCyclicPrefix ifftSignal];

Nbin = nobins * (NFFT + CPL);				        ## Total number of OFDM symbols

channelInput = reshape(SignalWithCyclicPrefix.', 1, Nbin);

signalLength = length(channelInput);

######### RICIAN ###############
Kdb = 2;  %Ratio of direct path and scattered paths
Mi = 1;   %Interpolation factor

K = 10 ^ (Kdb/10);
const = 1/(2*(K+1));
x = randn(1,signalLength);
y = randn(1,signalLength);
r = sqrt(const*((x + sqrt(2*K)).^2 + y.^2));
rt = zeros(1,Mi*length(r));
ki = 1;

for b=1:length(r)
  rt(ki:b*Mi) = r(b);
  ki = ki+Mi;
end
h = rt;

##################################

errors = [];

for noise = Noise
  
  awgnnoise = sqrt(noise/2)*randn(1, Nbin) + i*sqrt(noise/2)*randn(1, Nbin);

  channel = sqrt((NFFT+CPL)/NFFT)*h.*channelInput+ awgnnoise;
  channel = channel./h;
  
  receivedSignal = reshape(channel.', NFFT+CPL, nobins).';
  receivedSignal = receivedSignal(:, end-NFFT+1:end);
  receivedSymbols = (sqrt(DSC)/NFFT)*fftshift(fft(receivedSignal.')).';
  

  data = qamdemod(receivedSymbols, modulation);  
  finalData = de2bi(data.', symbolbits, "left-msb");
  data = reshape(data.', 1, NFFT*nobins);

  errors = [errors nnz(inputSymbols-data')];

endfor

errors = errors/N;

semilogy(SNRdb,errors);

legend("OFDM Simulation");

