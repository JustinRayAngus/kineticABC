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
Ecc = hdf5read(thisFile,'Ecc');
Ece = hdf5read(thisFile,'Ece');
f0 = hdf5read(thisFile,'f0');
Te0 = hdf5read(thisFile,'Te0');
tns = hdf5read(thisFile,'tout');
%display(Ecc);
%display(Ece);
display(Te0);

close(figure(1)); f1=figure(1);
plot(Ecc,sqrt(Ecc).*f0(:,1)); xlabel('\epsilon [eV]'); ylabel('\epsilon^1^/^2f_0 [1/eV]');
title('Initial EEDF');
hold on; plot(Ecc,sqrt(Ecc).*f0(:,2),'rx');

mom0 = sum(sqrt(Ecc).*f0(:,1))*(Ecc(2)-Ecc(1)); % should be one
mom2 = sum(Ecc.*sqrt(Ecc).*f0(:,1))*(Ecc(2)-Ecc(1)); % should be 3/2*Te

% filePath = '../build/';
% fileName = 'outputTest.h5';
% thisFile = [filePath,fileName];
% fileinfo = hdf5info(thisFile);
% var = hdf5read(thisFile,'ExtendibleArray')';
% display(var);