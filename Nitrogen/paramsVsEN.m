%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot some parameters for 2Term sims vs E/N
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

thisMVpm = 1:1:10;


for k=1:length(thisMVpm)

%%%   load my simulation results
%
filePath = ['./',num2str(thisMVpm(k)),'MVpm/'];
fileName = 'output.h5';
thisFile = [filePath,fileName];
fileinfo = hdf5info(thisFile);
Ecc = hdf5read(thisFile,'Ecc');
Ece = hdf5read(thisFile,'Ece');
F0 = hdf5read(thisFile,'F0');
zeroMom = hdf5read(thisFile,'zeroMom');
Te0 = hdf5read(thisFile,'Te');
t  = hdf5read(thisFile,'tout');
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
%display(Ecc);
%display(Ece);
%display(Te);
nt = length(F0(1,:));


Qelm = hdf5read(thisFile,'Qelm'); % [m^2]
Qmom = hdf5read(thisFile,'Qmom');
Qexc = hdf5read(thisFile,'Qexc');
Uexc = hdf5read(thisFile,'Uexc');



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
Qmom = Qelm+sum(Qexc,2)+sum(Qizn,2);

F0ce(2:length(Ece)) = 10.^interp1(log10(Ecc),log10(F0(:,nt)),log10(Ece(2:length(Ece))),'pchirp');
F0ce(1) = F0(1,nt);
Qmomcc = 10.^interp1(log10(Ece(2:length(Ece))),log10(Qmom(2:length(Ece))),log10(Ecc),'pchirp');
dF0dE = zeros(size(Ecc));
for i = 1:length(Ecc)
    deltaE = Ece(i+1)-Ece(i);
    dF0dE(i) = (F0ce(i+1)-F0ce(i))/deltaE;
    thisIntegrand(i) = -gamma/3.0*Ecc(i)^1.5*dF0dE(i) ...
                     /(sqrt(Ecc(i))*Qmomcc(i)+nunet(nt)/(Ng*gamma));
 %   thisIntegrand(i) = -gamma/3.0*Ecc(i)*dF0dE(i)/Qmomcc(i);
end

Te(k) = Te0(nt);
EN(k) = Ez/Ng*1e21;  % reduced E [Td]
muN(k) = sum(thisIntegrand*deltaE); % reduced mobility [1/s/m/V]
Vdrift(k) = muN(k)*EN(k)*1e-21*100;       % drift speed [cm/s]
alpha(k) = nunet(nt)/Vdrift(k);        % townsend coefficient [1/cm]
alphaN(k) = alpha(k)/Ng*1e6;           % reduced townsend coefficient [cm^2]
VT(k)     = 4.19e7*sqrt(Te(k));    % thermal speed [cm/s]
display(EN(k));
%display(alphaN);
%display(Vdrift);



end


%%%   plot normalized variables vs E/N
%
close(figure(1));
f1 = figure(1); set(f1,'position',[300 900 1500 420]);
%
subplot(1,3,1);
plot(EN,1.5*Te);
xlabel('reduced electric field [Td]');
ylabel('3/2 T_e [eV]'); grid on;
title('mean electron energy');
%
subplot(1,3,2);
semilogy(EN,alphaN);
xlabel('reduced electric field [Td]');
ylabel('reduced townsend coefficient [cm^2]');
title('reduced townsend coefficient');
grid on; grid off; grid on;
axis([0 500 1e-24 1e-16]);
%
subplot(1,3,3);
semilogy(EN,alphaN.*Vdrift);
xlabel('reduced electric field [Td]');
ylabel('rate constant [cm^3/s]');
title('ionization rate constant');
grid on; grid off; grid on;
axis([0 500 1e-18 1e-8]);



%%%   plot variables vs electric field
%
close(figure(2));
f2 = figure(2); set(f2,'position',[300 400 1500 420]);
%
subplot(1,3,1);
plot(thisMVpm,1.5*Te);
xlabel('electric field [MV/m]');
ylabel('3/2 T_e [eV]'); grid on;
title('mean electron energy at P=760 Torr');
%
subplot(1,3,2);
semilogy(thisMVpm,alpha);
xlabel('electric field [MV/m]');
ylabel('townsend coefficient [1/cm]');
title('townsend coefficient at P=760 Torr');
grid on; grid off; grid on;
%set(gca,'YMinorGrid','on')
axis([0 10 1e-2 1e3]);
set(gca,'xtick',0:2:10);
set(gca,'ytick',10.^(-2:1:3));
%
subplot(1,3,3);
semilogy(thisMVpm,alpha.*Vdrift);
xlabel('electric field [MV/m]');
ylabel('rate [1/s]');
title('ionization rate at P=760 Torr');
grid on; grid off; grid on;
%set(gca,'YMinorGrid','on')
axis([0 10 1e4 1e11]);
set(gca,'ytick',10.^(4:1:11));


