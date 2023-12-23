% Run this script for generating simulated CaV and BK open probabilities of
% the 1:1 BK-CaV complex in response to a voltage step from -80 mV to 0 mV,
% obtained from the 6-state ODE model (Eqs. S6-S11), the simplified
% Hodgkin-Huxley-type model (Eq. 23 coupled with Eq. 20)
% and the corresponding model assuming instantaneous activation (Eq. 29
% with k=1, m_{CaV}=m_{CaV,\infty}, coupled with Eq. 20), as shown in Figure 1BC. 
% Also, compare the results with 
% those obtained from the 70-state Markov chain model (Cox, 2014), 
% and the 6-state Markov chain model.

clear all 
clc
close all

% define global parameters
global k1_0 K1  k2_0 K2 n1 n2 alpha0 beta0 V_0  ca_c  n_c  k_minus k_plus N_tot beta1 alpha1 beta1_0 alpha1_0 rho

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
V=0; % voltage value  [mV]
V_0=-80; %initial voltage value [mv]

n_c=1; % number of calcium channels (CaVs) coupled with BK
% ca_o= -6/20*V+19;

N_tot=1; % the probabilities sum is 1 (N_tot=cx+cy+ox+oy+bx+by)

% window time at which it is applied the step voltage from -80 mV to V 
tv=[0 100];
time_in=-10;
time_f=100;

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
% ca_o= -6/20*V+19;

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

k1=k1_0*exp(-alpha0*V); % w^{+} defined in the paper 

k2=k2_0*exp(-beta0*V); % w^{-} defined in the paper 

fco_Ca_c=(ca_c^n1)/(ca_c^n1+K1^n1); % f^{+}(Ca_c) defined in the paper 
fco_Ca_o=(ca_o^n1)/(ca_o^n1+K1^n1); % f^{+}(Ca_o) defined in the paper 

foc_Ca_c=(K2^n2)/(K2^n2+ca_c^n2); % f^{-}(Ca_c) defined in the paper 
foc_Ca_o=(K2^n2)/(K2^n2+ca_o^n2); % f^{-}(Ca_o) defined in the paper 


k1_ca_o=k1*fco_Ca_o; % k^{+}_{o} defined in the paper 
k1_ca_c=k1*fco_Ca_c; % k^{+}_{c} defined in the paper 
k2_ca_o=k2*foc_Ca_o; % k^{-}_{o} defined in the paper 
k2_ca_c=k2*foc_Ca_c; % k^{-}_{c} defined in the paper 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BK-CaV model with n_c=1 non-inactivated CaV%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan=time_in:1e-2:time_f;
x_0=zeros(1, 2*(n_c+1)) ;
x_0(1)=1;

% 4-state ODE model for the BK activation
[T, Y]=ode15s('bk_n_cav_states_model_scale',tspan, x_0,[],tv,V_in);
bk_y=Y(:,n_c+2:end);
y_sum=sum(bk_y,2); % BK open prob, p_y, with non-inactivated CaV

figure (3)
hold on 
grid on
plot (T,Y)
plot(T,y_sum,'-b','linewidth',2);
legend('cx','ox','cy','oy','p_y')
xlim([0 20])
title('4-state ODE model for the BK activation in 1:1 BK-CaV complex with non-inactivated CaV')
xlabel('time [ms]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BK-CaV  model approximation with n_c=1 non-inactivated CaV%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yk1_0=0;
yk2_0=[0,1];

% 1-ODE model for the BK activation (Eq. 29 with k=1) with m_CaV= m_{CaV,\infty} 
[T1, y1qssminf]=ode45('bk_cav_1state_model_scale_qss_minf',tspan, yk1_0,[],tv,V_in);

% 1-ODE model for the BK activation (Eq. 23) and 1-ODE model for the CaV activation 
[T1, Y2]=ode45('bk_cav_2state_model_scale_qss',tspan, yk2_0,[],tv,V_in);
yk2=Y2(:,1); % m_{BK} 
c2=Y2(:,2);  % c in the paper for CaV

figure (4)
grid on
hold on
plot(T,y_sum,'-b') 
plot(T1,yk2,'--r')
plot(T1,y1qssminf,'-.g')

% m_BK at ss for V value
m_BK_inf=k1_ca_o.*m_inf.*(alpha+beta_scale+k2_ca_c)./((k1_ca_o+k2_ca_o).*(alpha+k2_ca_c)+beta_scale.*k2_ca_c);

tf=20; % [ms]
idx_tf=find(T1==tf);
plot(T1(idx_tf),m_BK_inf,'o')
xlim([0 tf])
xlabel('time [ms]');
title('1:1 BK-CaV complex with non-inactivated CaV')
legend('p_y','m_{BK} with m_{CaV}', 'm_{BK_{approx.}}','m_{BK,\infty}')

%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CaV model %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
tspan1=time_in:1e-2:time_f;
x_10=[1 0];
k_minus=2*10^(-3); % gamma in the paper
k_plus=2.5*10^(-3); % delta/Ca_{CaV} in the paper

% 2-state ODE model for CaV
[T2c, Y2c]=ode45('cav_2states_model_scale',tspan1, x_10,[],tv,V_in);
c=Y2c(:,1);
o=Y2c(:,2);
b=1-c-o;
h_ca=1-b;

% 1-state ODE model for CaV
[T1, b1]=ode45('cav_1state_h_scale',tspan1, 0,[],tv,V_in);

h1=1-b1; % h for CaV with 1-state ODE model
 
figure (5)
hold on 
grid on
plot(T2c,c,'-c',T2c,h_ca,'--g',T2c,b,'--r',T2c,o,'--b')
plot(T1,h1,'--k')
legend('c (2-ODE)','h (2-ODE)','b (2-ODE)','o (2-ODE)','h (1-ODE)')
title('CaV model')
xlim([0 100])
xlabel('time [ms]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BK-CaV model with n_c=1 inactivating CaV %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan6= time_in:1e-2:time_f;
x_0=[0, 0,0,0,0 ];
% 6-state ODE model for 1:1 BK-CaV
[T6, Y6]=ode45('bk_cav_6states_model_scale',tspan6, x_0,[],tv,V_in);

oy=Y6(:,1);
cy=Y6(:,2);
ox=Y6(:,3);
bx=Y6(:,4);
by=Y6(:,5);
cx=N_tot-oy-cy-ox-bx-by; 

y=oy+cy+by; % p_y, BK open prob with the 6-state ODE model
figure (6)
hold on 
grid on
plot (T6,oy,'b',T6,cy,'r',T6,ox,'g',T6,cx,'m',T6,bx,'k',T6,by,'c');
legend( 'oy','cy','ox','cx','bx','by')
title ('6-state ODE model for 1:1 BK-CaV')
xlim([0 100])
xlabel('time [ms]');

% load the simulation obtained from the Cox model
load('MC_sim/Cox_model.mat')

% load the simulation obtained from the 3*2-state Markov chain model 
load('MC_sim/6_state_model.mat')

f=figure (1);
hold on 
grid on
y_a_tot_1=y1qssminf.*h1; % BK open prob with 1-state ODE for m_{BK} (Eq. 29 with k=1), m_{CaV}=m_{CaV,\infty}, and 1-state ODE model for h 
y_a_tot_2=yk2.*h1; % BK open prob with 1-state ODE for m_{BK} (Eq. 23), 1-state ODE for m_{CaV}, and 1-state ODE model for h 
plot(t_cox,mean(open_BK_cox_model),'-.','color',[0.5,0.5,0.5]); % Cox model results
plot(t_6_state-0.01,mean(open_BK_6_state_model),'-.k') % 6-state MC results, time shift, in 6-state MC simulations voltage step is applied t=0.01

plot(T6,y,'-b','linewidth',1); % 
plot(T,y_a_tot_2,'--r','linewidth',1);
plot(T,y_a_tot_1,'-.g','linewidth',1);
ylim([0 0.4])
xlim([0 20])
xlabel('time [ms]');
ylabel('BK open prob.')
legend('MC model by Cox', '6-state MC model', '6-state ODE model', 'm_{BK}h', 'm_{BK_{approx.}}h (m_{CaV}=m_{CaV,\infty})')
title('Figure 1C')
%set(f,'Position',[10 10 150 125])


f=figure(2);
hold on
grid on
plot(t_cox,mean(open_CaV_cox_model),'-.','color',[0.5,0.5,0.5]); % Cox model results
plot(t_6_state-0.01,mean(open_CaV_3_state_model),'--k') % 3-state MC results
plot(T,o,'-b','linewidth',1);
plot(T,m_inf*h1,'-.g','linewidth',1);
xlim([0 20])
ylim([0 0.66])
legend('MC model by Cox', '3-state MC model', '3-state ODE model', 'm_{CaV,\infty}h (1-ODE for h)')
ylabel('CaV open prob.')
xlabel('time [ms]')
title('Figure 1B')
%set(f,'Position',[10 10 150 125])


