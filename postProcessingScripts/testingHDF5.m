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
data = hdf5read(thisFile,'Ece');
display(data);