%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   Helium
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


%%%   load my simulation results
%
filePath = './EN10/';
%filePath = '../test2/';
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
Qizn = hdf5read(thisFile,'Qizn');
Qelm = hdf5read(thisFile,'Qelm');
Qexc = hdf5read(thisFile,'Qexc');
Uizn = hdf5read(thisFile,'Uizn');
nunet = hdf5read(thisFile,'nunet');
W = hdf5read(thisFile,'W'); % W = W(:,2);
D = hdf5read(thisFile,'D'); % D = D(:,2);
Ez = hdf5read(thisFile,'E');      % [V/m]
mM = hdf5read(thisFile,'mM');
Ng = hdf5read(thisFile,'Ng');     % [1/m^3]
EN = Ez/Ng*1e21;  % reduced E [Td]
%display(Ecc);
%display(Ece);
%display(Te);
nt = length(F0(1,:));

%%%   plot EEDF
%
close(figure(3)); f3=figure(3); set(f3,'position',[0 100 1000 400]);
subplot(1,2,1);
semilogy(Ecc,F0(:,1));
xlabel('\epsilon [eV]'); ylabel('F_0 [1/eV^3^/^2]');
title(['He4 EEDF evolution for E/N = ',num2str(EN,3),' Td']);
hold on; plot(Ecc,F0(:,round(nt/12)),'magenta');
hold on; plot(Ecc,F0(:,round(nt/2)),'color',[0 0.5 0]);
hold on; plot(Ecc,F0(:,nt),'r');
%hold on; plot(E0,F00,'black--'); % plot steady state solution
axis([0 max(Ece) 1e-18 1]);
%formatSpec = '%10.2e\n';
formatSpec = 2;
lg1=legend('t=0',['t=',num2str(1e9*t(round(nt/12)),2),' ns'], ...
       ['t=',num2str(1e9*t(round(nt/2)),formatSpec),' ns'], ...
       ['t=',num2str(1e9*t(nt),formatSpec),' ns']); 
set(lg1','location','best');

mom0 = sum(sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be one
mom2 = sum(Ecc.*sqrt(Ecc).*F0(:,1))*(Ecc(2)-Ecc(1)); % should be 3/2*Te

dFlux = zeros(length(Ecc),length(Flux(1,:)));
deltaE = Ecc(2)-Ecc(1);
for i = 1:length(Ecc)
    dFlux(i,:) = (Flux(i+1,:)-Flux(i,:))/(deltaE*sqrt(Ecc(i)));
end
close(figure(112));
figure(112); %hold on; 
semilogx(Ecc,dFlux(:,nt),'b',Ecc,ExcS(:,nt)+IznS(:,nt)-nunet(nt)*F0(:,nt),'r--');

legend('div Flux','Source'); 
xlabel('\epsilon [eV]');

Qelm = hdf5read(thisFile,'Qelm'); % [m^2]
Qmom = hdf5read(thisFile,'Qmom');
Qexc = hdf5read(thisFile,'Qexc');
Uexc = hdf5read(thisFile,'Uexc');


figure(3);
subplot(1,2,2);
tns = 1e9*t;
plot(tns,zeroMom,'b',tns,Te,'r'); 
xlabel('t [ns]'); ylabel('Moments');
axis([0 max(tns) 0 1.1*max(Te)]);
legend('zero', 'Te [eV]');
title('Evolution of moments');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%   post calculate transport coefficients
%%%
%%%

%%%   compute mobility 
%
econst = 1.6022e-19;
meconst = 9.1094e-31;
gamma = sqrt(2*econst/meconst);
Qmom = Qelm+Qexc+Qizn;

F0ce(2:length(Ece)) = 10.^interp1(log10(Ecc),log10(F0(:,nt)),log10(Ece(2:length(Ece))),'pchirp');
F0ce(1) = F0(1,nt);
Qmomcc = 10.^interp1(log10(Ece(2:length(Ece))),log10(Qmom(2:length(Ece))),log10(Ecc),'pchirp');
dF0dE = zeros(size(Ecc));
for i = 1:length(Ecc)
    dF0dE(i) = (F0ce(i+1)-F0ce(i))/deltaE;
    thisIntegrand(i) = -gamma/3.0*Ecc(i)^1.5*dF0dE(i) ...
                     /(sqrt(Ecc(i))*Qmomcc(i)+nunet(nt)/(Ng*gamma));
 %   thisIntegrand(i) = -gamma/3.0*Ecc(i)*dF0dE(i)/Qmomcc(i);
end
muN = sum(thisIntegrand*deltaE); % reduced mobility [1/s/m/V]
Vdrift = muN*EN*1e-21*100;       % drift speed [cm/s]
alpha = nunet(nt)/Vdrift;        % townsend coefficient [1/cm]
alphaN = alpha/Ng*1e6;           % reduced townsend coefficient [cm^2]
VT     = 4.19e7*sqrt(Te(nt));    % thermal speed [cm/s]
display(alphaN);
display(Vdrift);