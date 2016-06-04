%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%      script to run comsol with matlab function
%%%
%%%      NOTE: must open this script after using command prompt
%%%      to open COMSOl with MATLAB by first going to ../bin/maci64
%%%      and then typing './comsol server matlab'
%%% 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   NOTE that when opening matlab as desribed above, the default
%%%   plot values are changed. So I manually put them back to what I
%%%   want (see startup.m file)

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

set(0,'defaultlinelinewidth',1)
load('cmap2.mat')
set(0,'DefaultFigureColormap',white_to_red_cmap);

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

already_loaded = 1;
if(already_loaded==0)
%    clear all;
    %%%       set path and filename to COMSOL program
    %
    filepath = '../test2/';
    filename = 'N2_test2';  % COMSOL file
    model = mphload([filepath,filename]);
end

%%%       obtain solution
%
%model.sol('sol1').run;


%%%       get stuff to plot
%
info = mphsolinfo(model, 'soltag','sol1');
EN = info.solvals;              % reduced electric field [Td]
Te = mphglobal(model,'be.Te'); 

Emax=str2num(model.param.get('Emax'));
E = 0:0.05:Emax; % xdomain to interpret COMSOL results to

close(figure(1)); figure(1);
F0 = mphinterp(model,{'f'},'coord',E,'Dataset','dset1')'; % [1/eV^1.5]
semilogy(E,F0); xlabel('\epsilon [eV]'); ylabel('F_0 [1/eV^3^/^2]');
axis([0 30 1e-10 1]); title('EEDF N_2')
legend(['E/N = ',num2str(EN,3),' Td']);

N2comsol.Te = Te;
N2comsol.EN = EN;
N2comsol.E = E;
N2comsol.F0 = F0;
save('../test2/N2comsolResults.mat','N2comsol');

