filename = "errors.txt";

SNRdb = [1:10];

try
  errorList = [];
  legendNameList = [""];
  load(filename);
catch
  errorList = [];
  legendNameList = [""];
end_try_catch



figure(1);

hold on;

for i = 1:3
  semilogy(SNRdb, errorList(i,:), legendNameList(i,:));
endfor

hold off;



figure(2);

hold on;

for i = 4:6
  semilogy(SNRdb, errorList(i,:), legendNameList(i,:));
endfor

hold off;



figure(3);

hold on;

for i = 7:7
  semilogy(SNRdb, errorList(i,:), legendNameList(i,:));
endfor

hold off;




figure(4);

hold on;

for i = 8:12
  semilogy(SNRdb, errorList(i,:), legendNameList(i,:));
endfor


hold off;
