% Run this script for generating simulated CaV and BK open probabilities
% obtained from the 3*n*2-state Markov chain model, as shown in Figure 1BC (n=1) and
% Figure 2C (n=4), where n (n_c in the code) is the number of CaV channels 
% in 1:n BK-CaV complex

close all
clear all
V=0; %voltage value
delta_t=.01; % time step
n_c=1; % number of calcium channels coupled with BK for Figure 1BC
% n_c=4; % number of calcium channels coupled with BK for Figure 2C


% cav kinetics
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

% Ca2+ concentration at the internal mouth of the CaV (Ca_{CaV}) when the
% CaV is open
ca_o_rca=i_ca/(8*pi*D_ca*F*conv_F*r_ca)*exp(-r_ca/(sqrt(D_ca/k_B/B_tot)))*conv_microM;

% Ca2+ concentration at the BK when the CaV is open
ca_o=i_ca/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;

k_minus=2*10^(-3); % gamma in the paper
k_plus=2.5*10^(-3); % delta/Ca_{CaV} in the paper

ca_c=0.2; % background Ca2+ concentration
ca_d=ca_o_rca;

%CaV parameters
cv_par=[1.2979    1.0665   -0.0639    0.0703    0.3090]; % dark

alpha1_0=cv_par(1);
beta1_0=cv_par(2);
alpha1=cv_par(3);
beta1=cv_par(4);
rho=cv_par(5);

alpha=alpha1_0*exp(-alpha1*V); % alpha in the paper
beta_s=beta1_0*exp(-beta1*V);
beta=(beta_s+alpha)* rho; % beta in the paper

% transition probabilities for CaV
p_cav_co=alpha*delta_t;
p_cav_oc=beta*delta_t;  

p_cav_ob=k_plus*ca_d*delta_t;
p_cav_bo=k_minus*delta_t;

p_cav_cc=1-p_cav_co;
p_cav_oo=1-p_cav_oc - p_cav_ob;
p_cav_bb=1-p_cav_bo;


% BK parameters
bk_par=[  1.1093e3   3.3206e3  2.3298    0.0223    1.6012   16.5793    0.1000    0.4614];

k1_0=bk_par(1)/1000;

k2_0=bk_par(2)/1000;

K1=bk_par(6);

K2=bk_par(7);

n1=bk_par(3);

n2=bk_par(8);

beta0=bk_par(4);

alpha0=-bk_par(5)*beta0;

k1=k1_0*exp(-alpha0*V); % w^{+} defined in the paper

k2=k2_0*exp(-beta0*V); % w^{-} defined in the paper

n_sim=1000; % number of the simulations
tot_dt=10000; % time step number
cav_c=zeros(n_c,n_sim,tot_dt); % c
cav_o=zeros(n_c,n_sim,tot_dt); % o
cav_b=zeros(n_c,n_sim,tot_dt); % b
bk=zeros(n_sim,tot_dt); % p_y (open probability for the BK)
Ca_v=zeros(n_c,n_sim,tot_dt); % Ca2+ concentration at the BK for 1 CaV 


% Monte Carlo simulations
for num_sim=1:n_sim
    c=1;o=0;b=0;Ca=0.2;
    cav_c(:, num_sim,1)=c;
    cav_o(:,num_sim,1)=o;
    cav_b(:,num_sim,1)=b;
    bk(num_sim,1)=0;
    Ca_v(:,num_sim,1)=Ca;
    x0=0;
    for num_it=2:tot_dt
        
        for num_ch=1:n_c
            y=rand(1,1);
            
            c=cav_c(num_ch, num_sim,num_it-1);
            o=cav_o(num_ch,num_sim,num_it-1);
            b= cav_b(num_ch,num_sim,num_it-1);
            
            if c==1
                
                if  y<p_cav_cc
                    c=1;o=0;b=0;
                    Ca=0;
                else
                    c=0;o=1;b=0;
                    Ca=ca_o;
                end
                
                
            elseif o==1
                if  y<p_cav_oo
                    c=0;o=1;b=0;
                    Ca=ca_o;
                elseif y>p_cav_oo && y<1-p_cav_oc
                    c=0;o=0;b=1;
                    Ca=0;
                else
                    c=1;o=0;b=0;
                    Ca=0;
                end
                
                
            elseif b==1
                
                if  y<p_cav_bb
                    c=0;o=0;b=1;
                    Ca=0;
                else
                    c=0;o=1;b=0;
                    Ca=ca_o;
                end
            end
            
        
            
            cav_c(num_ch,num_sim,num_it)=c;
            cav_o(num_ch,num_sim,num_it)=o;
            cav_b(num_ch,num_sim,num_it)=b;
            
            Ca_v(num_ch,num_sim,num_it)=Ca;
        
        end
        Ca_n_c=sum(Ca_v(:,num_sim,num_it));
        if Ca_n_c==0
            Ca_n_c=ca_c;
        end
        % transition probabilities for BK
        fco_Ca=k1*(Ca_n_c^n1)/(Ca_n_c^n1+K1^n1)*delta_t;
        foc_Ca=k2*(K2^n2)/(K2^n2+Ca_n_c^n2)*delta_t;
        
        fcc_Ca=1-fco_Ca;
        foo_Ca=1-foc_Ca;
        
        
        x=rand(1,1);
        
        if x0==0
            if x<fcc_Ca
                bk(num_sim,num_it)=0;
            else
                bk(num_sim,num_it)=1;
            end
        else
            if x<foc_Ca
                bk(num_sim,num_it)=0;
            else
                bk(num_sim,num_it)=1;
            end
        end
        
        x0=bk(num_sim,num_it);
    end

end

t_sim=0.01:0.01:100;

p_bk_mean=mean(bk);
figure
grid on
hold on
plot(t_sim,p_bk_mean,'k')
xlabel('time [ms]');
ylabel('BK open prob.')



cav1_o=squeeze(cav_o(1,:,:));
cav1_b=squeeze(cav_b(1,:,:));
cav1_c=squeeze(cav_c(1,:,:));
Ca_v1=squeeze(Ca_v(1,:,:));

p_cav1_c_mean=mean(cav1_c);
p_cav1_o_mean=mean(cav1_o);
p_cav1_b_mean=mean(cav1_b);
p_Ca_v1_mean=mean(Ca_v1);

% cav2_o=squeeze(cav_o(2,:,:));
% cav3_o=squeeze(cav_o(3,:,:));
% cav4_o=squeeze(cav_o(4,:,:));
% 
% cav2_b=squeeze(cav_b(2,:,:));
% cav3_b=squeeze(cav_b(3,:,:));
% cav4_b=squeeze(cav_b(4,:,:));
% 
% cav2_c=squeeze(cav_c(2,:,:));
% cav3_c=squeeze(cav_c(3,:,:));
% cav4_c=squeeze(cav_c(4,:,:));
% 
% Ca_v2=squeeze(Ca_v(2,:,:));
% Ca_v3=squeeze(Ca_v(3,:,:));
% Ca_v4=squeeze(Ca_v(4,:,:));
% 
% p_Ca_v2_mean=mean(Ca_v2);
% p_Ca_v3_mean=mean(Ca_v3);
% p_Ca_v4_mean=mean(Ca_v4);
% 
% p_cav2_c_mean=mean(cav2_c);
% p_cav2_o_mean=mean(cav2_o);
% p_cav2_b_mean=mean(cav2_b);
% 
% p_cav3_c_mean=mean(cav3_c);
% p_cav3_o_mean=mean(cav3_o);
% p_cav3_b_mean=mean(cav3_b);
% 
% p_cav4_c_mean=mean(cav4_c);
% p_cav4_o_mean=mean(cav4_o);
% p_cav4_b_mean=mean(cav4_b);

figure
grid on, hold on
plot(t_sim,p_cav1_o_mean,'k')
legend('o')
% plot(t_sim,p_cav1_o_mean,'r',t_sim,p_cav2_o_mean,'b',t_sim,p_cav3_o_mean,'g',t_sim,p_cav4_o_mean,'m')
% legend('o1','o2','o3','o4')
xlabel('time [ms]');
ylabel('CaV open prob.')

% figure
% grid on, hold on
% plot(t_sim,p_cav1_c_mean,'c')
% legend('c')
% xlabel('time [ms]');
% plot(t_sim,p_cav1_c_mean,'r',t_sim,p_cav2_c_mean,'b',t_sim,p_cav3_c_mean,'g',t_sim,p_cav4_c_mean,'m')
% legend('c1','c2','c3','c4')

% 
% % 
% figure
% grid on, hold on
% plot(t_sim,p_cav1_b_mean,'r')
% legend('b')
% xlabel('time [ms]');
% plot(t_sim,p_cav1_b_mean,'r',t_sim,p_cav2_b_mean,'b',t_sim,p_cav3_b_mean,'g',t_sim,p_cav4_b_mean,'m')
% legend('b1','b2','b3','b4')

% 
% figure
% grid on, hold on
% plot(t_sim,p_Ca_v1_mean,'r')
% legend('Ca')
% xlabel('time [ms]');
% plot(t_sim,p_Ca_v1_mean,'r',t_sim,p_Ca_v2_mean,'b',t_sim,p_Ca_v3_mean,'g',t_sim,p_Ca_v4_mean,'m')
% legend('Ca1','Ca2','Ca3','Ca4')




