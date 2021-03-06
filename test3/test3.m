%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   compare N2 results with steady state values from COMSOL
%%%   for E/N = 1000 Td 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

sharing = 'Equal';  % energy sharing for ionization procesess
%sharing = 'Zero';

%%%   load comsol results 
%
load('comsolResultsN2_test3.mat');
E0  = comsolResults.E;
EN0 = comsolResults.EN;

if(strcmp(sharing,'Equal'))
    Te0 = comsolResults.Te_equal;
    F00 = comsolResults.F0equal;
    fileName = 'outputEqualSharing.h5';
elseif(strcmp(sharing,'Zero'))
    Te0 = comsolResults.Te_zero;
    F00 = comsolResults.F0zero;
    fileName = 'outputZeroSharing.h5';
else
    display('sharing is either Equal or Zero!!!');    
end

filePath = './';
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
close(figure(3)); fg=figure(3); set(fg,'position',[0 100 1000 400]);
subplot(1,2,1);
semilogy(Ecc,F0(:,1));
xlabel('\epsilon [eV]'); ylabel('F_0 [1/eV^3^/^2]');
title(['N2 EEDF evolution for E/N = ',num2str(EN,3),' Td and ', sharing,' energy sharing']);
hold on; plot(Ecc,F0(:,round(nt/12)),'magenta');
hold on; plot(Ecc,F0(:,round(nt/8)),'color',[0 0.5 0]);
hold on; plot(Ecc,F0(:,nt),'r');
hold on; plot(E0,F00,'black--'); % plot steady state solution
axis([0 100 1e-6 0.1]);
%formatSpec = '%10.2e\n';
formatSpec = 2;
lg1=legend('t=0',['t=',num2str(1e9*t(round(nt/12)),2),' ns'], ...
       ['t=',num2str(1e9*t(round(nt/8)),formatSpec),' ns'], ...
       ['t=',num2str(1e9*t(nt),formatSpec),' ns'],'Soln'); 
set(lg1','location','best');




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

dFlux = zeros(length(Ecc),length(Flux(1,:)));
deltaE = Ecc(2)-Ecc(1);
for i = 1:length(Ecc)
    dFlux(i,:) = (Flux(i+1,:)-Flux(i,:))/(deltaE*sqrt(Ecc(i)));
end
close(figure(111));
figure(111); semilogx(Ecc,dFlux(:,nt),'b',Ecc,ExcS(:,nt)+IznS(:,nt)-nunet(nt)*F0(:,nt),'r--');
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
hold on; plot(tns,0*tns+Te0,'black--');
xlabel('t [ns]'); ylabel('Moments');
axis([0 max(tns) 0 1.1*max(Te)]);
legend('zero', 'Te [eV]','Te Soln','location','best');
title('Evolution of moments');