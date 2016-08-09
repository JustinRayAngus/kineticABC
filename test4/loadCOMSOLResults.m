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
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

set(0,'defaultlinelinewidth',1)
load('cmap2.mat')
set(0,'DefaultFigureColormap',white_to_red_cmap);

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = './';
filename = './N2_test3';
model = mphload([filepath,filename]);


info = mphsolinfo(model, 'soltag','sol1');
EN = info.solvals;
Emax = str2double(model.param.get('Emax'));
% nE = str2double(model.param.get('nE'));
% dE = Emax/nE;
x = 0:Emax/2000:Emax;

%%%   get the EEDF
%
F0equal = mphinterp(model,{'f'},'coord',x,'Dataset','dset1')'; % [eV]
F0zero  = mphinterp(model,{'f'},'coord',x,'Dataset','dset2')';
Te_equal = mphglobal(model,'be.Te','Dataset','dset1');
Te_zero = mphglobal(model,'be.Te','Dataset','dset2');

close(figure(1)); figure(1);
semilogy(x,F0equal,'b',x,F0zero,'r');
legend('equal','zero');

comsolResults.EN = EN;
comsolResults.Te_equal = Te_equal;
comsolResults.Te_zero = Te_zero;
comsolResults.F0equal = F0equal;
comsolResults.F0zero = F0zero;
comsolResults.E = x;

save('comsolResultsN2_test3.mat','comsolResults');

