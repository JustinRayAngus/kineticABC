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
F0 = hdf5read(thisFile,'F0');
Te0 = hdf5read(thisFile,'Te0');
tns = hdf5read(thisFile,'tout');
Qelm = hdf5read(thisFile,'Qelm'); % [m^2]
Ng = hdf5read(thisFile,'Ng'); % 1/m^3
Flux = hdf5read(thisFile,'Flux');
W = hdf5read(thisFile,'W');
D = hdf5read(thisFile,'D');
%display(Ecc);
%display(Ece);
%display(Te0);

close(figure(2)); f2=figure(2);
plot(Ecc,sqrt(Ecc).*F0(:,1)); xlabel('\epsilon [eV]'); ylabel('\epsilon^1^/^2F_0 [1/eV]');
title('Initial EEDF');
hold on; plot(Ecc,sqrt(Ecc).*F0(:,101),'r--');

mom0 = sum(sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be one
mom2 = sum(Ecc.*sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be 3/2*Te

close(figure(3)); 
f3=figure(3);
plot(Ece,Flux(:,2),'black'); 
hold on; plot(Ece,Flux(:,90),'r--');
%hold on; plot(Ece,W(:,2).*sqrt(Ece),'r');
%hold on; plot(Ece,D(:,2),'b');
xlabel('\epsilon [eV]'); ylabel('Flux');
title('RHS Flux term');

% dFlux = zeros(size(Ecc));
% for i = 1:length(Ecc)
%     dFlux(i) = Flux(i+1,2)-Flux(i,2);
% end
% 
% close(figure(4)); 
% f4=figure(4);
% plot(Ecc,-dFlux/(Ecc(2)-Ecc(1)),'r--'); xlabel('\epsilon [eV]'); 
% ylabel('-dFlux');
% title('RHS Flux term');

% filePath = '../build/';
% fileName = 'outputTest.h5';
% thisFile = [filePath,fileName];
% fileinfo = hdf5info(thisFile);
% var = hdf5read(thisFile,'ExtendibleArray')';
% display(var);