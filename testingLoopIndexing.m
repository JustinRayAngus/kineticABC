%%%
%%%
%%%
%%%
clear all;

U = 5.2;
Emax = 10;
nE = 10;
dE = Emax/10;
Ecc = dE/2:dE:Emax-dE/2;
Ece = [0 Ecc+dE/2];

%%% step 1 is to find index for cell that contains U
%
[~,kU] = min(abs(Ecc-U));


%%% step 2 is to loop over E, and find index for cell that contains
%%% E+U and 2E+U
%
thisEU  = zeros(size(Ecc)); %  Ecc + U
this2EU = zeros(size(Ecc)); % 2Ecc + U
kEU  = zeros(size(Ecc));    % index for cell that contains  Ecc+U
k2EU = zeros(size(Ecc));    % index for cell that contains 2Ecc+U 


%%% first get stuff for i=1
%
thisEU(1) = Ecc(1)+U;
this2EU(1) = 2*Ecc(1)+U;
if(thisEU(1)<Ece(kU+1))
    kEU(1)  = kU;
else
    kEU(1) = kU+1;
end
k2EU(1) = kU+1;


%%% 2nd loop through all other cells
%
for i=2:length(Ecc)
    thisEU(i) = Ecc(i)+U;
    this2EU(i) = 2*Ecc(i)+U;
    kEU(i)  = kEU(i-1)+1;
    k2EU(i) = k2EU(i-1)+2;
    
    if(this2EU(i)>=max(Ece))
        this2EU(i) = max(Ece);
        k2EU(i) = nE;
    end
    if(thisEU(i)>=max(Ece))
        thisEU(i) = max(Ece);
        kEU(i) = nE;
        break;
    end
    
end
    
        
        
        
        