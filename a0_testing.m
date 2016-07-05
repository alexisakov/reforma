%% Start IRIS if it is not already running:
try, irisversion;, catch, addpath 'C:\Users\AIsakov\Downloads\IRIS'; irisstartup;, end;

%% CLEAR: clean the workspace
clear all; close all;
%%
clear all;
m1 = model('reforma.mod','linear',true);
m1=solve(m1);
m1 = sstate(m1,'growth=',true,'blocks=',true);
chksstate(m1)
get(m1,'sstate')