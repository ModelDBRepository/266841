%Run this script for simulating the CaV and BK open probabilities  in
%response to a voltage step (V in the code) (or the simulated action potential (load('data_AP.mat'))), 
%obtained from the 70-state Markov chain model (Cox, 2014) of the 1:1
%BK-CaV complex.

clear all
close all
clc

tsim = 100; %ms
dt = 0.01;  %ms
g_ca = 2.8; %pS
Vca = 60;   % mV
Vstep = 0;  % mV (voltage step value)
kpiuB = 500; % microM-1*s-1
Ca_bg = 0.1; % microM
Btotal = 30; % microM
F = 96485;
Dca = 250;   % microm*s-1;
rep1 = 0;
rep2 = 0;
t_end=rep1+rep2+tsim;

% voltage step and time
 V=[ones((rep1)/dt,1)*-80; ones((tsim+dt)/dt,1)*Vstep;  ones((rep2)/dt,1)*-80];
 Voltage_time=0:dt:t_end;

% %%% voltage step and time
% t1=-0.5:0.01:0;
% t2=0:0.01:20;
% t3=20.01:0.01:24;
% t_end=t3(end);
% V=[ones(length(t2),1)*Vstep;  ones(length(t3),1)*-80];
% Voltage_time=0:dt:t_end;

%%% load the simulated action potential
% load('data_AP.mat')
% Voltage_time=time_AP;
% V=v_AP;
% t_end=Voltage_time(end);

tot_dt=length(Voltage_time);
n_sim=1000;
Yopen_CaV=zeros(n_sim,tot_dt);
open_stoc= zeros(n_sim,tot_dt);
closed_stoc= zeros(n_sim,tot_dt);
 for k=1:n_sim
 %% Stochastic Ca  
 X0=[1 0 0 0 0 0 0];
 [Y,t] = MonteCarlo_seven_states(Voltage_time,V,X0,t_end,struct,dt);
 Yopen=Y(6,:);
 Yopen_CaV(k,:)=Yopen;
 % Ca concentration
 r = 0.013; % microM;
 ica=(g_ca*(V-Vca))*1e-3; %pA
 Caf = (abs(ica)*1e9)/(8*pi*Dca*F*r)*exp(-r/sqrt(Dca/(kpiuB*Btotal)))+Ca_bg; 
 struct.Ca=Caf.*Yopen'; 
 struct.Voltage_time=[0:dt:t_end]';
 struct.Voltage=V;
  
%% stochastic BK
X0=[0 0 0 0 0 1 0 0 0 0];
[X,t] = BK_stochastic(X0,t_end,dt,struct);
open_stoc(k,:)=X(1,:)+X(2,:)+X(3,:)+X(4,:)+X(5,:);
closed_stoc(k,:)=X(6,:)+X(7,:)+X(8,:)+X(9,:)+X(10,:);
end

figure
plot(t,Caf)
title('Caloc')

% plot stochastic opening of Ca channel
figure
plot(t,mean(Yopen_CaV,1))
title('Open stochastic CaV')

% plot stochastic opening of BK channel
figure
plot(t,mean(open_stoc,1))
title('Open stochastic BK')

figure
plot(t,mean(closed_stoc,1)+mean(open_stoc,1))
title('sum stochastic')

