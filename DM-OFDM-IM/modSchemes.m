
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







