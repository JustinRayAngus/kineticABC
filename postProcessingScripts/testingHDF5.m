%%%
%
%   looking at hdf5 output file
%
%%%
clear all;

filePath = '../build/';
fileName = 'output.h5';
thisFile = [filePath,fileName];
fileinfo = hdf5info(thisFile);
Ecc = hdf5read(thisFile,'Ecc')';
Ece = hdf5read(thisFile,'Ece')';
f0 = hdf5read(thisFile,'f0')';
%display(Ecc);
%display(Ece);

close(figure(1)); f1=figure(1);
plot(Ecc,sqrt(Ecc).*f0); xlabel('\epsilon [eV]'); ylabel('\epsilon^1^/^2f_0 [1/eV]');
title('Initial EEDF');

mom0 = sum(sqrt(Ecc).*f0)*(Ecc(2)-Ecc(1)); % should be one
mom2 = sum(Ecc.*sqrt(Ecc).*f0)*(Ecc(2)-Ecc(1)); % should be 3/2*Te