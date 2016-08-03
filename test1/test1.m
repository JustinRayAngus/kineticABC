%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   test1 computes EEDF using elastic momentum tranfer only, which has
%%%   a semi-analytical solution
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

filePath = './';
fileName = 'output.h5';
thisFile = [filePath,fileName];
fileinfo = hdf5info(thisFile);
Ecc = hdf5read(thisFile,'Ecc');
Ece = hdf5read(thisFile,'Ece');
F0 = hdf5read(thisFile,'F0');
zeroMom = hdf5read(thisFile,'zeroMom');
Te = hdf5read(thisFile,'Te');
t   = hdf5read(thisFile,'tout');
Flux = hdf5read(thisFile,'Flux');
ExcS = hdf5read(thisFile,'ExcS');
IznS = hdf5read(thisFile,'IznS');
W = hdf5read(thisFile,'W'); 
D = hdf5read(thisFile,'D');


nt = length(F0(1,:));
close(figure(3)); f3=figure(3); set(f3,'position',[0 100 1000 400]);
subplot(1,2,1); 
semilogy(Ecc,F0(:,1)); xlabel('\epsilon [eV]'); ylabel('F_0 [1/eV^3^/^2]');
title('EEDF evolution');
hold on; plot(Ecc,F0(:,round(nt/4)),'magenta');
hold on; plot(Ecc,F0(:,round(nt/2)),'color',[0 0.5 0]);
hold on; plot(Ecc,F0(:,nt),'r');

mom0 = sum(sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be one
mom2 = sum(Ecc.*sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be 3/2*Te


Ez = hdf5read(thisFile,'E');      % [V/m]
mM = hdf5read(thisFile,'mM');
Ng = hdf5read(thisFile,'Ng');     % [1/m^3]
Qelm = hdf5read(thisFile,'Qelm'); % [m^2]
Qmom = hdf5read(thisFile,'Qmom');

exponent = zeros(size(Qelm));
deltaE = Ece(2)-Ece(1);
for i = 2:length(Qelm)
    thisval = -6*Ng^2*mM/Ez^2*Ece(i)*Qelm(i)^2*deltaE;
    exponent(i) = exponent(i-1)+thisval;
end
FSoln = exp(exponent);
normC = sum(sqrt(Ece).*FSoln*deltaE);
FSoln = FSoln/normC;

figure(3);
%VE4 = Ez^2/(3*Ng^2*max(Qelm)^2)/mM;
%VE2 = sqrt(VE4);
%A = 25.6909/(2*pi*VE2).^1.5; 
%FSoln = A*exp(-(Ecc/VE2).^2);
%sum(sqrt(Ecc).*Fsoln.*(Ecc(2)-Ecc(1)))
TeSoln = 2/3*sum(sqrt(Ece.^3).*FSoln.*deltaE);
figure(3); hold on; plot(Ece,FSoln,'black--');
formatSpec = '%10.2e\n';
legend('t=0',['t=',num2str(t(round(nt/4)),formatSpec)], ...
       ['t=',num2str(t(round(nt/2)),formatSpec)], ...
       ['t=',num2str(t(nt),formatSpec)],'Soln'); 
   %axis([0 60 0 1.1*max(sqrt(Ece).*FSoln)]);

figure(3);
subplot(1,2,2);
plot(t,zeroMom,'b',t,Te,'r'); 
hold on ;plot(t,0*t+TeSoln,'black--');
xlabel('t [s]'); ylabel('Moments');
axis([0 max(t) 0 1.1*max(Te)]);
legend('zero', 'Te', 'Soln');
title('Evolution of moments');
