clear; clc; 
%% Time trial for Algorithm 11.2 
close all;
time112 = tic;
tic;
NR_Burgers_CNRKW3_FD;
toc; 
disp(time112);

%% Time trial for Algorithm 11.3
close all;
time113 = tic;
tic;
NR_Burgers_CNRKW3_FD_RS;
toc;
disp(time113);

%% Time comparison
ratio = time112 / time113;
disp(ratio);
save("timetrial.mat", "time112", "time113", "ratio");