% Run this script for generating simulated CaV and BK open probabilities of
% the 1:n BK-CaV complex (n=4) in response to a voltage step from -80 mV to 0 mV,
% obtained from the 2*(n+1)ODE model for the BK activation (Eqs. S19-S24) coupled with 1-state
% ODE model for h (Eq. 20), 
% from the reduced ODE model considering CaV activation kinetics for the BK 
% activation (Eq. 26) coupled with 1-state ODE model for h (Eq. 20), 
% and from the simplification assuming m_{CaV}=m_{CaV,\infty} for the BK 
% activation (Eq. 29) coupled with 1-state ODE model for h (Eq. 20), as shown in Figure 2C.
% Also, compare the results with those obtained from the 3*n*2-state Markov 
% chain model. 

clear all 
clc
close all


% define global parameters
global k1_0 K1  k2_0 K2 n1 n2 alpha0 beta0 V_0 ca_c N_tot n_c  k_minus k_plus  beta1 alpha1 beta1_0 alpha1_0 rho

% BK parameters
bk_par=[  1.1093   3.3206  2.3298    0.0223    1.6012   16.5793    0.1000    0.4614];

k1_0=bk_par(1);

k2_0=bk_par(2);

K1=bk_par(6);

K2=bk_par(7);

n1=bk_par(3);

n2=bk_par(8);

beta0=bk_par(4);

alpha0=-bk_par(5)*beta0;

ca_c=0.2; % background Ca2+ concentration 
V=0; % voltage value
V_0=-80; %initial voltage value

n_c=4; % number of calcium channels coupled with BK
% ca_o= -6/20*V+19;
N_tot=1;

% window time at which it is applied the step voltage from -80 mV to V 
tv=[0 100];

% voltage step from -80 mV to V value and then back to -80 mV (at t=100 ms)
V_in=[V -80]; 


%%% compute calcium level at 7nm of the pore (cav sensor) and at 13nm (bk sensor)
D_ca=250; %microm^2 s^-1
F=9.6485*10^4; %C mol^-1
conv_F=10^(-15);  % mol a M/microm^3
conv_microM=10^6; % M to microM
r_bk=13*10^(-3);% microm
r_ca=7*10^(-3); %microm
k_B=500; %microM^-1 * s^-1
B_tot=30; %microM
Eca=60; %mV
conv_V=10^(-3); %mV to V
g_ca=2.8; %pS 
conv_S=10^(-12); %pS to S
i_ca=abs((V-Eca))*conv_V*g_ca*conv_S; % C/sec

ca_o_rca=i_ca/(8*pi*D_ca*F*conv_F*r_ca)*exp(-r_ca/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_o=i_ca/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;

% CaV parameters
cv_par=[1.2979    1.0665   -0.0639    0.0703    0.3090]; 
alpha1_0=cv_par(1);

beta1_0=cv_par(2);

alpha1=cv_par(3);

beta1=cv_par(4);

rho=cv_par(5);

alpha=alpha1_0*exp(-alpha1*V); % alpha in the paper
beta=beta1_0*exp(-beta1*V);
beta_scale=(beta+alpha)* rho; % beta in the paper

m_inf=alpha/(alpha+beta_scale); % m_{CaV,\infty} 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2*(n+1) ODE model for 1:n BK-CaV complex with n non-inactivated CaVs (Eq. S19-S24) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan= -10:1e-2:100;
x_0=zeros(1, 2*(n_c+1)) ;
x_0(1)=1;
[T, Y]=ode45('bk_n_cav_states_model_scale',tspan, x_0,[],tv,V_in);
bk_y=Y(:,n_c+2:end);
y4_sum=sum(bk_y,2); % BK open prob with 4 non-inactivated CaVs (n=4)

figure
hold on 
grid on
plot (T,Y)
plot(T,y4_sum,'-r','linewidth',2);
legend('ccccx','cccox','ccoox','cooox','oooox','ccccy','cccoy','ccooy','coooy','ooooy','p_y')
xlim([0 20])
xlabel('time [ms]');
title('2*(4+1) ODE model for 1:4 BK-CaV complex with non-inactivated CaVs')

n_c=1;
x_0=zeros(1, 2*(n_c+1)) ;
x_0(1)=1;
[T, Y]=ode45('bk_n_cav_states_model_scale',tspan, x_0,[],tv,V_in);
bk_y=Y(:,n_c+2:end);
y1_sum=sum(bk_y,2); % BK open prob with 1 non-inactivated CaV 

n_c=2;
x_0=zeros(1, 2*(n_c+1)) ;
x_0(1)=1;
[T, Y]=ode45('bk_n_cav_states_model_scale',tspan, x_0,[],tv,V_in);
bk_y=Y(:,n_c+2:end);
y2_sum=sum(bk_y,2); % BK open prob with 2 non-inactivated CaVs 

n_c=3;
x_0=zeros(1, 2*(n_c+1)) ;
x_0(1)=1;
[T, Y]=ode45('bk_n_cav_states_model_scale',tspan, x_0,[],tv,V_in);
bk_y=Y(:,n_c+2:end);
y3_sum=sum(bk_y,2); % BK open prob with 3 non-inactivated CaVs 

n_c=4;
%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CaV model %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
tspan1=-10:1e-2:100;
x_10=[1 0];
k_minus=2*10^(-3);
k_plus=2.5*10^(-3);

% 2-state ODE model for CaV
[T1c, Y1c]=ode15s('cav_2states_model_scale',tspan1, x_10,[],tv,V_in);
c=Y1c(:,1);
o=Y1c(:,2);
b=1-c-o;
%h_ca=1-b;

% 1-state ODE model for CaV
[T1, b1]=ode45('cav_1state_h_scale',tspan1, 0,[],tv,V_in);
h_ca=1-b1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1:n BK-CaV model approximation with n=1, 2, 3, 4 non-inactivated CaVs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximate m_BK (Eq. 29) and assume m_{CaV}=m_{CaV,\infty} 
yk1_0=0;
[T1, y1minf]=ode45('bk_cav_1state_model_scale_qss_minf',tspan, yk1_0,[],tv,V_in);
[T1, y2minf]=ode45('bk_n2_cav_1state_model_scale_approx_qss_minf',tspan, yk1_0,[],tv,V_in);
[T1, y3minf]=ode45('bk_n3_cav_1state_model_scale_approx_qss_minf',tspan, yk1_0,[],tv,V_in);
[T1, y4minf]=ode45('bk_n4_cav_1state_model_scale_approx_qss_minf',tspan, yk1_0,[],tv,V_in);

% Assume CaV activation kinetics  (Eq. 26 for m_{BK})
yk1_0=[0 0];
[T1, z1]=ode45('bk_cav_2state_model_scale_qss',tspan, yk1_0,[],tv,V_in);
[T1, z2]=ode45('bk_n2_cav_2state_model_scale_approx_qss',tspan, yk1_0,[],tv,V_in);
[T1, z3]=ode45('bk_n3_cav_2state_model_scale_approx_qss',tspan, yk1_0,[],tv,V_in);
[T1, z4]=ode45('bk_n4_cav_2state_model_scale_approx_qss',tspan, yk1_0,[],tv,V_in);

yk1=z1(:,1);
yk2=z2(:,1);
yk3=z3(:,1);
yk4=z4(:,1);


figure
hold on
grid on
lw=1;
plot(T,yk1,'--c','linewidth',lw)
plot(T, yk2,'--g','linewidth',lw)
plot(T,yk3,'--m','linewidth',lw)
plot(T, yk4,'--r','linewidth',lw)
plot(T,y1minf,'-.c','linewidth',lw)
plot(T,y2minf,'-.g','linewidth',lw)
plot(T,y3minf,'-.m','linewidth',lw)
plot(T,y4minf,'-.r','linewidth',lw)
plot(T,y1_sum,'-c','linewidth',lw)
plot(T,y2_sum,'-g','linewidth',lw)
plot(T,y3_sum,'-m','linewidth',lw)
plot(T,y4_sum,'-r','linewidth',lw)
title('BK activation in complexes with n non-inactivated CaVs (n=1, 2, 3, 4)'); 
xlim([0 20])
xlabel('time [ms]');
legend('m_{BK}^{(1)}',' m_{BK}^{(2)}',' m_{BK}^{(3)}',' m_{BK}^{(4)}', 'm_{BK_{approx.}}^{(1)} (m_{CaV,\infty})','m_{BK_{approx.}}^{(2)} (m_{CaV,\infty})',' m_{BK_{approx.}}^{(3)} (m_{CaV,\infty})',' m_{BK_{approx.}}^{(4)} (m_{CaV,\infty})', 'm_{BK}^{(1)} (full ODE)', 'm_{BK}^{(2)} (full ODE)', 'm_{BK}^{(3)} (full ODE)', 'm_{BK}^{(4)} (full ODE)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1:n BK-CaV complex with n inactivating CaVs (n=4) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_1=nchoosek(n_c,1)*h_ca.*(1-h_ca).^(n_c-1);
h_2=nchoosek(n_c,2)*h_ca.^2.*(1-h_ca).^(n_c-2);
h_3=nchoosek(n_c,3)*h_ca.^3.*(1-h_ca).^(n_c-3);
h_4=h_ca.^n_c;
yk1_h=yk1.*h_1;
yk2_h=yk2.*h_2;
yk3_h=yk3.*h_3;
yk4_h=yk4.*h_4;

y1m_h=y1minf.*h_1;
y2m_h=y2minf.*h_2;
y3m_h=y3minf.*h_3;
y4m_h=y4minf.*h_4;

y1s_h=y1_sum.*h_1;
y2s_h=y2_sum.*h_2;
y3s_h=y3_sum.*h_3;
y4s_h=y4_sum.*h_4;

figure
hold on
grid on
plot(T,h_ca,'--g',T,h_1,'--m',T,h_2,'--c',T,h_3,'--r',T,h_4,'--b')
legend('h','4h(1-h)^3','6h^2(1-h)^2','h^3(1-h)','h^4')
xlabel('time [ms]');

% load the simulation obtained from the 3*4*2-state Markov chain model 
load('MC_sim/3_4n_2_state_model.mat')


f=figure;
hold on 
grid on
y_s_tot=y1s_h+y2s_h+y3s_h+y4s_h;
y_a_tot=yk1_h+yk2_h+yk3_h+yk4_h;
y_b_tot=y1m_h+y2m_h+y3m_h+y4m_h;
plot(t_sim-0.01,mean(open_BK_3_n4_2_state_model),'-.k') % 3*4*2-state Markov chain model results - time shift, in 3*4*2-state MC simulations voltage step is applied t=0.01
plot(T,y_s_tot,'-b','linewidth',1);
plot(T,y_a_tot,'--r','linewidth',1);
plot(T,y_b_tot,'-.g','linewidth',1);
%set(f,'position',[20 20 150 125])
xlim([0 20])
ylim([0 0.8])
xlabel('time [ms]');
ylabel('BK open prob.')
legend ('3*n*2-state MC model', '2*(n+1)-state ODE model for BK activation','m_{BK}^{(n)}','m_{BK_{approx}}^{(n)} with m_{CaV,\infty}')
title('Figure 2C  (n=4)')


