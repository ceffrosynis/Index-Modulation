pkg load communications;
clc; clear all;

modulation = 16;

symbols = qammod(0:(modulation-1), modulation);

figure(1);

title("mod-1");

indexSymbolA = real(symbols) + 3 == 0 | real(symbols) + 1 == 0;
indexSymbolB = real(symbols) - 3 == 0 | real(symbols) - 1 == 0;
hold on;
scatter(real(symbols(indexSymbolA)), imag(symbols(indexSymbolA)));
scatter(real(symbols(indexSymbolB)), imag(symbols(indexSymbolB)));
hold off;

figure(2);

title("mod-2");

indexSymbolA = real(symbols) + 3 == 0 | real(symbols) - 1 == 0;
indexSymbolB = real(symbols) - 3 == 0 | real(symbols) + 1 == 0;
hold on;
scatter(real(symbols(indexSymbolA)), imag(symbols(indexSymbolA)));
scatter(real(symbols(indexSymbolB)), imag(symbols(indexSymbolB)));
hold off;

figure(3);

title("mod-3");

indexSymbolA = imag(symbols) > 0 & (real(symbols) + 3 == 0 | real(symbols) - 1 == 0) | imag(symbols) < 0 & (real(symbols) + 1 == 0 | real(symbols) - 3 ==0);
indexSymbolB = imag(symbols) < 0 & (real(symbols) + 3 == 0 | real(symbols) - 1 == 0) | imag(symbols) > 0 & (real(symbols) + 1 == 0 | real(symbols) - 3 ==0);
hold on;
scatter(real(symbols(indexSymbolA)), imag(symbols(indexSymbolA)));
scatter(real(symbols(indexSymbolB)), imag(symbols(indexSymbolB)));
hold off;

figure(4);

title("mod-4");

indexSymbolA = (imag(symbols) == -1 | imag(symbols) == 3) & (real(symbols) == -3 | real(symbols) == 1) | (imag(symbols) == -3 | imag(symbols) == 1) & (real(symbols) == -1 | real(symbols) == 3);
indexSymbolB = ~indexSymbolA;
hold on;
scatter(real(symbols(indexSymbolA)), imag(symbols(indexSymbolA)));
scatter(real(symbols(indexSymbolB)), imag(symbols(indexSymbolB)));
hold off;


figure(5);

title("mod-5");

indexSymbolA = (real(symbols) == -3 | real(symbols) == 1) & (imag(symbols) == 3 | imag(symbols) == -1);
indexSymbolB = (real(symbols) == -1 | real(symbols) == 3) & (imag(symbols) == 3 | imag(symbols) == -1);
indexSymbolC = (real(symbols) == -3 | real(symbols) == 1) & (imag(symbols) == 1 | imag(symbols) == -3);
indexSymbolD = (real(symbols) == -1 | real(symbols) == 3) & (imag(symbols) == 1 | imag(symbols) == -3);
hold on;
scatter(real(symbols(indexSymbolA)), imag(symbols(indexSymbolA)));
scatter(real(symbols(indexSymbolB)), imag(symbols(indexSymbolB)));
scatter(real(symbols(indexSymbolC)), imag(symbols(indexSymbolC)));
scatter(real(symbols(indexSymbolD)), imag(symbols(indexSymbolD)));
hold off;


figure(6);

title("mod-6");

groupedSymbols = generateSymbolGroups (16);
hold on;
for i = 1:size(groupedSymbols, 1);
  scatter(real(groupedSymbols(i,:)), imag(groupedSymbols(i,:)), 100, "filled");
endfor
legend("mode 1", "mode 2", "mode 3", "mode 4", "location", "east");
hold off;

figure(7);

title("mod-7");

groupedSymbols = generateSymbolGroups (32);
hold on;
for i = 1:size(groupedSymbols, 1);
  scatter(real(groupedSymbols(i,:)), imag(groupedSymbols(i,:)), 100, "filled");
endfor
legend("mode 1", "mode 2", "mode 3", "mode 4", "mode 5", "mode 6", "mode 7", "mode 8", "location", "east");
zoom out;
zoom out;
hold off;

figure(8);

title("mod-8");

groupedSymbols = generateSymbolGroups (64);
hold on;
for i = 1:size(groupedSymbols, 1);
  scatter(real(groupedSymbols(i,:)), imag(groupedSymbols(i,:)), 100, "filled");
endfor
legend("mode 1", "mode 2", "mode 3", "mode 4", "mode 5", "mode 6", "mode 7", "mode 8", "mode 9", "mode 10", "mode 11", "mode 12", "mode 13", "mode 14", "mode 15", "mode 16", "location", "east");
hold off;






