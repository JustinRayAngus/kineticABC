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
zeroMom = hdf5read(thisFile,'zeroMom');
Te0 = hdf5read(thisFile,'Te0');
t   = hdf5read(thisFile,'tout');
Flux = hdf5read(thisFile,'Flux');
W = hdf5read(thisFile,'W');
D = hdf5read(thisFile,'D');
%display(Ecc);
%display(Ece);
%display(Te0);

nt = length(F0(1,:));
close(figure(2)); f2=figure(2); set(f2,'position',[0 100 1000 400]);
subplot(1,2,1); 
plot(Ecc,sqrt(Ecc).*F0(:,1)); xlabel('\epsilon [eV]'); ylabel('\epsilon^1^/^2F_0 [1/eV]');
title('EEDF evolution');
hold on; plot(Ecc,sqrt(Ecc).*F0(:,round(nt/8)),'magenta');
hold on; plot(Ecc,sqrt(Ecc).*F0(:,round(nt/4)),'color',[0 0.5 0]);
hold on; plot(Ecc,sqrt(Ecc).*F0(:,nt),'r');

mom0 = sum(sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be one
mom2 = sum(Ecc.*sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be 3/2*Te

% close(figure(3)); 
% f3=figure(3);
% plot(Ece,Flux(:,2),'black'); 
% hold on; plot(Ece,Flux(:,nt),'r--');
% %hold on; plot(Ece,W(:,2).*sqrt(Ece),'r');
% %hold on; plot(Ece,D(:,2),'b');
% xlabel('\epsilon [eV]'); ylabel('Flux');
% title('RHS Flux term');

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

Ez = hdf5read(thisFile,'E');      % [V/m]
mM = hdf5read(thisFile,'mM');
Ng = hdf5read(thisFile,'Ng');     % [1/m^3]
Qelm = hdf5read(thisFile,'Qelm'); % [m^2]
Qexc = hdf5read(thisFile,'Qexc');

exponent = zeros(size(Qelm));
deltaE = Ece(2)-Ece(1);
for i = 2:length(Qelm)
    thisval = -6*Ng^2*mM/Ez^2*Ece(i)*Qelm(i)^2*deltaE;
    exponent(i) = exponent(i-1)+thisval;
end
FSoln = exp(exponent);
normC = sum(sqrt(Ece).*FSoln*deltaE);
FSoln = FSoln/normC;

%VE4 = Ez^2/(3*Ng^2*max(Qelm)^2)/mM;
%VE2 = sqrt(VE4);
%A = 25.6909/(2*pi*VE2).^1.5; 
%FSoln = A*exp(-(Ecc/VE2).^2);
%sum(sqrt(Ecc).*Fsoln.*(Ecc(2)-Ecc(1)))
TeSoln = 2/3*sum(sqrt(Ece.^3).*FSoln.*deltaE);
figure(2); hold on; plot(Ece,sqrt(Ece).*FSoln,'black--');
formatSpec = '%10.2e\n';
legend('t=0',['t=',num2str(t(round(nt/8)),formatSpec)], ...
       ['t=',num2str(t(round(nt/4)),formatSpec)], ...
       ['t=',num2str(t(nt),formatSpec)],'Soln'); axis([0 60 0 1.1*max(sqrt(Ece).*FSoln)]);

figure(2);
subplot(1,2,2);
plot(t,zeroMom,'b',t,Te0,'r'); 
hold on ;plot(t,0*t+TeSoln,'black--');
xlabel('t [s]'); ylabel('Moments');
axis([0 max(t) 0 1.1*max(Te0)]);
legend('zero', 'Te', 'Soln');
title('Evolution of moments');
