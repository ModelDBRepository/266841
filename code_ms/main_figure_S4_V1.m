% Run this script for simulating the response of  the 1:1 BK-CaV complex to
% different voltage steps. In particular, simulate the responses
% obtained from the 6-state ODE model (Eqs. S6-S11), the simplified Hodgkin-Huxley-type
% model (Eq. 23 coupled with Eq. 20)
% and the further simplified model assuming instantaneous activation  (Eq. 29
% with k=1, m_{CaV}=m_{CaV,\infty}, coupled with Eq. 20),
% as shown in Figure S4. Also, compare the results with 
% those obtained from the 70-state Markov chain model (Cox, 2014).


clear all 
clc
close all

plot_cox=1; % if 1, plot the results obtained from the Cox model 

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
V_0=-80; %initial voltage value

n_c=1; % number of calcium channels (CaVs) coupled with BK
% ca_o= -6/20*V+19;

n_cav_tot=1000;
n_bk_tot=1000;
conv_na=1e-12*1e-3*1e9;
conv_pa=1e-12*1e-3*1e12;
N_tot=1; % the probabilities sum is 1 (N_tot=cx+cy+ox+oy+bx+by)

% CaV parameters
cv_par=[1.2979    1.0665   -0.0639    0.0703    0.3090]; % 

alpha1_0=cv_par(1);

beta1_0=cv_par(2);

alpha1=cv_par(3);

beta1=cv_par(4);

rho=cv_par(5);

V_v=[-40 -20 0 20 40]; % step voltage values

for V=V_v
    
    % window time at which it is applied the step voltage from -80 mV to V
    tv=[0 20];
    time_in=-10;
    time_f=24;
    
    % voltage step from -80 mV to V value and then back to -80 mV (at t=20 ms)
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
   % E_K=-75; %mV
    E_K=-96;  
   conv_V=10^(-3); %mV to V
    g_ca=2.8; %pS
    g_bk=100;
    conv_S=10^(-12); %pS to S
    i_ca=abs((V-Eca))*conv_V*g_ca*conv_S; % C/sec
    
    ca_o_rca=i_ca/(8*pi*D_ca*F*conv_F*r_ca)*exp(-r_ca/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
    ca_o=i_ca/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
    % ca_o= -6/20*V+19;
    
    % CaV parameters
    alpha=alpha1_0*exp(-alpha1*V); % alpha in the paper
    beta=beta1_0*exp(-beta1*V);
    beta_scale=(beta+alpha)* rho; % beta in the paper  
    m_inf=alpha/(alpha+beta_scale); % m_{CaV,\infty}
    
    alphaV0=alpha1_0*exp(-alpha1*V_0); 
    betaV0=beta1_0*exp(-beta1*V_0);
    beta_scaleV0=(betaV0+alphaV0)* rho;
    m_inf_V0=alphaV0/(alphaV0+beta_scaleV0); % m_{CaV,\infty}
    
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
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BK-CaV  model approximation with n_c=1 non-inactivated CaV%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yk1_0=0;
    yk2_0=[0,1];
    
    % 1-ODE model for the BK activation (Eq. 29 with k=1) with m_CaV= m_{CaV,\infty}
    [T1, y1qssminf]=ode45('bk_cav_1state_model_scale_qss_minf',tspan, yk1_0,[],tv,V_in);
    
    % 1-ODE model for the BK activation (Eq. 23) and 1-ODE model for the CaV activation
    [T1a, Y2]=ode45('bk_cav_2state_model_scale_qss',tspan, yk2_0,[],tv,V_in);
    
    yk2=Y2(:,1); % m_{BK}
    c2=Y2(:,2);  % c in the paper for CaV
   
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
     
    y_a_tot_1=y1qssminf.*h1; % BK open prob with 1-state ODE for m_{BK} (Eq. 29 with k=1), m_{CaV}=m_{CaV,\infty}, and 1-state ODE model for h
    y_a_tot_2=yk2.*h1; % BK open prob with 1-state ODE for m_{BK} (Eq.23), 1-state ODE for m_{CaV}, and 1-state ODE model for h
    
    % step voltage
    V_vect=V*ones(length(T),1); 
    idx_T=find(T==20);
    idx_T0=find(T==0);
    V_vect(idx_T+1:end)=V_0;
    V_vect(1:idx_T0-1)=V_0;
    
    wl=150;
    wh=150;
    T_in=-1;
    minf_h_CaV=[m_inf_V0*h1(1:idx_T0-1); m_inf*h1(idx_T0:idx_T); m_inf_V0*h1(idx_T+1:end)];
    
    f=figure(1);
    hold on
    grid on
    plot(T,V_vect,'-b','linewidth',1);
    xlim([T_in 24])
    set(f,'Position',[10 10 wl wh])
    ylim([-82 42])
    title('Figure S4A')
    xlabel('time [ms]');
    ylabel('Voltage steps [mV]')
    
    f=figure(2);
    hold on
    grid on
    plot(T, conv_pa*n_cav_tot*g_ca*o.*(V_vect-Eca),'-b','linewidth',1);
    xlim([T_in 24])
    set(f,'Position',[10 10 wl wh])
    ylim([-186 0])
    xlabel('time [ms]');
    ylabel('I_{CaV} [pA]')  
    title('Figure S4C')

    f=figure(3);
    hold on
    grid on
    plot(T,conv_pa*n_cav_tot*g_ca*minf_h_CaV.*(V_vect-Eca),'-g','linewidth',1);
    xlim([T_in 24])
    set(f,'Position',[10 10 wl wh])
    ylim([-186 0]) 
    xlabel('time [ms]');
    ylabel('I_{CaV} [pA]')  
    title('Figure S4D')

    f=figure (4);
    hold on
    grid on
    plot(T6,conv_na*n_bk_tot*g_bk*y.*(V_vect-E_K),'-b','linewidth',1); %
    xlim([T_in 24])
    ylim([0 7])
    set(f,'Position',[10 10 wl wh])
    xlabel('time [ms]');
    ylabel('I_{BK} [nA]')
    title('Figure S4F')
    
    f=figure (5);
    hold on
    grid on
    plot(T,conv_na*n_bk_tot*g_bk*y_a_tot_2.*(V_vect-E_K),'-r','linewidth',1);
    xlim([T_in 24])
    ylim([0 80])
    ylim([0 7])
    set(f,'Position',[10 10 wl wh])
    xlabel('time [ms]');
    ylabel('I_{BK} [nA]')
    title('Figure S4G')
    
    f=figure(6);
    hold on
    grid on
    plot(T,conv_na*n_bk_tot*g_bk*y_a_tot_1.*(V_vect-E_K),'-g','linewidth',1);
    xlim([T_in 24])
    ylim([0 7])
    set(f,'Position',[10 10 wl wh])
    xlabel('time [ms]');
    ylabel('I_{BK} [nA]')
    title('Figure S4H')
    
end

if plot_cox
    
    t2=0:0.01:20;
    t3=20.01:0.01:24;
    t=[t2 t3];
    Vstep=0;
    V_v0=[ones(length(t2),1)*Vstep;  ones(length(t3),1)*-80];
    Vstep=-40;
    V_minus40=[ones(length(t2),1)*Vstep;  ones(length(t3),1)*-80];
    Vstep=-20;
    V_minus20=[ones(length(t2),1)*Vstep;  ones(length(t3),1)*-80];
    Vstep=20;
    V_20=[ones(length(t2),1)*Vstep;  ones(length(t3),1)*-80];
    Vstep=40;
    V_40=[ones(length(t2),1)*Vstep;  ones(length(t3),1)*-80];
    
    load ('MC_sim/Cox_model_voltage_steps.mat')
    f=figure;
    hold on
    grid on
    plot(t,conv_pa*n_cav_tot*g_ca*open_cav_cox_minus40.*(V_minus40-Eca)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_pa*n_cav_tot*g_ca*open_cav_cox_minus20.*(V_minus20-Eca)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_pa*n_cav_tot*g_ca*open_cav_cox_0.*(V_v0-Eca)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_pa*n_cav_tot*g_ca*open_cav_cox_20.*(V_20-Eca)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_pa*n_cav_tot*g_ca*open_cav_cox_40.*(V_40-Eca)',':','color',[0.5,0.5,0.5],'linewidth',1);
    ylim([-186 0])
    xlim([T_in 24])
    title('Figure S4B')
    xlabel('time [ms]');
    ylabel('I_{CaV} [pA]')  
    set(f,'Position',[10 10 wl wh])
    
    f=figure;
    hold on
    grid on
    plot(t,conv_na*n_bk_tot*g_bk*open_bk_cox_minus40.*(V_minus40-E_K)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_na*n_bk_tot*g_bk*open_bk_cox_minus20.*(V_minus20-E_K)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_na*n_bk_tot*g_bk*open_bk_cox_0.*(V_v0-E_K)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_na*n_bk_tot*g_bk*open_bk_cox_20.*(V_20-E_K)',':','color',[0.5,0.5,0.5],'linewidth',1);
    plot(t,conv_na*n_bk_tot*g_bk*open_bk_cox_40.*(V_40-E_K)',':','color',[0.5,0.5,0.5],'linewidth',1);
    ylim([0 7])
    xlim([T_in 24])
    set(f,'Position',[10 10 wl wh])
    title('Figure S4E')
    xlabel('time [ms]');
    ylabel('I_{BK} [nA]')
    
end
    
    
    
