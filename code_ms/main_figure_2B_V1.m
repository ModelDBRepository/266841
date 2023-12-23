% Run this script for generating steady-state BK activation functions and
% time-constants for BK channels in complexes with 1, 2 or 4 CaVs, as shown  
% in Figure 2B, where the BK activation is modeled by Eq. 26 and Eq. 29.

clear all 
clc
close all

% define global parameters
global k1_0 K1  k2_0 K2 n1 n2 alpha0 beta0 ca_o ca_o_rca ca_c 

% BK parameters
BK_par=[  1.1093   3.3206  2.3298    0.0223    1.6012   16.5793    0.1000    0.4614];

k1_0=BK_par(1);

k2_0=BK_par(2);

K1=BK_par(6);

K2=BK_par(7);

n1=BK_par(3);

n2=BK_par(8);

beta0=BK_par(4);

alpha0=-BK_par(5)*beta0;

ca_c=0.2; % background Ca2+ concentration 

%%% parameters for computing calcium level at 7nm of the pore (cav sensor) and at 13nm (bk sensor)
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

% CaV parameters
cv_par=[1.2979    1.0665   -0.0639    0.0703    0.3090]; %
alpha1_0=cv_par(1);

beta1_0=cv_par(2);

alpha1=cv_par(3);

beta1=cv_par(4);

rho=cv_par(5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute parameters by varying V  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V values
V_vect=-60:0.1:60;
count=0;
for idx_V=V_vect
     
	count=count+1;
    V=idx_V;
    
	%CaV parameters
    
    alpha=alpha1_0*exp(-alpha1*V);
    beta=beta1_0*exp(-beta1*V);
    beta_scale=(beta+alpha)*rho;
    m_inf_v(count)=alpha/(alpha+beta_scale);

	% computa Ca2+ concentration at CaV sensor and BK sensor 
    i_ca=abs((V-Eca))*conv_V*g_ca*conv_S; % C/sec
    ca_o_rca=i_ca/(8*pi*D_ca*F*conv_F*r_ca)*exp(-r_ca/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
    ca_o=i_ca/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
    %ca_o = -6/20*V+19;
	
    k1=k1_0*exp(-alpha0*V); % w^{+} 
    k2=k2_0*exp(-beta0*V);  % w^{-}
	
	% store CaV and BK parameters
    k1_v(count)=k1_0*exp(-alpha0*V); 
    k2_v(count)=k2_0*exp(-beta0*V);
    
    alpha_v(count)=alpha; 
    beta_v(count)=beta_scale;

    fco_Ca_c_v(count)=k1*(ca_c^n1)/(ca_c^n1+K1^n1); % k^{+}_{c}
    foc_Ca_c_v(count)=k2*(K2^n2)/(K2^n2+ca_c^n2); % k^{-}_{c}
    
    fco_Ca_o_v(count)=k1*(ca_o^n1)/(ca_o^n1+K1^n1); % k^{+}_{o}
    foc_Ca_o_v(count)=k2*(K2^n2)/(K2^n2+ca_o^n2);	% k^{-}_{o}	
    
    fco_2Ca_o_v(count)=k1*((2*ca_o)^n1)/((2*ca_o)^n1+K1^n1); % k^{+}_{o2}
    foc_2Ca_o_v(count)=k2*(K2^n2)/(K2^n2+(2*ca_o)^n2); % k^{-}_{o2}
    
    fco_3Ca_o_v(count)=k1*((3*ca_o)^n1)/((3*ca_o)^n1+K1^n1); % k^{+}_{o3}
    foc_3Ca_o_v(count)=k2*(K2^n2)/(K2^n2+(3*ca_o)^n2); % k^{-}_{o3}
    
    fco_4Ca_o_v(count)=k1*((4*ca_o)^n1)/((4*ca_o)^n1+K1^n1); % k^{+}_{o4}
    foc_4Ca_o_v(count)=k2*(K2^n2)/(K2^n2+(4*ca_o)^n2); % k^{-}_{o4}
    
end



k1_ca_o=fco_Ca_o_v;
k1_ca_c=fco_Ca_c_v;
k2_ca_o=foc_Ca_o_v;
k2_ca_c=foc_Ca_c_v;
k1_2ca_o=fco_2Ca_o_v;
k2_2ca_o=foc_2Ca_o_v;
k1_3ca_o=fco_3Ca_o_v;
k2_3ca_o=foc_3Ca_o_v;
k1_4ca_o=fco_4Ca_o_v;
k2_4ca_o=foc_4Ca_o_v;

m_inf=m_inf_v;
tau_v=1./(alpha_v+beta_v);
alpha=alpha_v;
beta=beta_v;

%%% n_CaV=1; 
% time constant (tau) and steady-state (ss) BK activation for equation (26) with
% k=n_ca=1 (for k=1 equation (26) becomes equation (23),
% where tau (tauBK_qss_1) and ss (mBKinf_qss_1) are given by equation (24)  
n_ca=1;

tauBK_qss_1=(alpha+beta+k2_ca_c)./((alpha+k2_ca_c).*(k1_ca_o+k2_ca_o)+beta.*k2_ca_c);
mBKinf_qss_1=k1_ca_o.*m_inf.*tauBK_qss_1;

% tau (tauBK_qss_minf_1) and ss (mBKinf_qss_minf_1) BK activation given by  
% equation (29) with k=n_ca=1;
ca1_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
ca1_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
tauBK_qss_minf_1=1./((1-m_inf).*k2_ca_c+m_inf_v.*(k1_ca_o+k2_ca_o));
mBKinf_qss_minf_1=k1_ca_o.*m_inf.*tauBK_qss_minf_1;

 
%%% n_CaV=2;
% tau (tauBK_qss_2) and ss (mBKinf_qss_2) BK activation for equation (26) 
% with k=n_ca=2 (see equation (S36)) 
n_ca=2;

A1=beta./(2*alpha+k2_ca_c+beta);
B_0=((alpha+k1_ca_o+k2_ca_o).*(1-A1)+2.*beta+k2_ca_c.*A1); % B_2 in the paper
A2=2*beta./B_0;
D2=2*(1-m_inf).*m_inf.*k1_ca_o./B_0;

tauBK_qss_2=1./(k2_ca_c.*A1.*A2+(k1_ca_o+k2_ca_o).*(1-A1).*A2 +(k1_2ca_o+k2_2ca_o).*(1-A2));
mBKinf_qss_2=(2*m_inf_v.*(1-m_inf).*k1_ca_o+k1_2ca_o.*m_inf.^n_ca - k2_ca_c.*D2.*A1-(k1_ca_o+k2_ca_o).*(1-A1).*(D2)+(k1_2ca_o+k2_2ca_o).*D2).*tauBK_qss_2;

% tau (tauBK_qss_minf_2) and ss (mBKinf_qss_minf_2) BK activation given by  
% equation (29) with k=n_ca=2;

n_ca=2;
ca2_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
ca2_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
ca2_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);

tauBK_qss_minf_2=1./((1-m_inf).^n_ca.*k2_ca_c+2*m_inf_v.*(1-m_inf).*(k1_ca_o+k2_ca_o)+m_inf_v.^n_ca.*(k1_2ca_o+k2_2ca_o));
mBKinf_qss_minf_2=( 2*m_inf_v.*(1-m_inf).*k1_ca_o+k1_2ca_o.*m_inf.^n_ca).*tauBK_qss_minf_2;


% %% n_CaV=4;
% tau (tauBK_qss_4) and ss (mBKinf_qss_4) BK activation for equation (26) 
% with k=n_ca=4 (see equation (S36)) 
n_ca=4;

A1=beta./(4*alpha+k2_ca_c+beta);
D1=0;
B2=((3*alpha+k1_ca_o+k2_ca_o).*(1-A1)+2*beta+k2_ca_c.*A1);
A2=2*beta./B2;
D2=4*m_inf.*(1-m_inf).^3.*k1_ca_o./B2; % k2_ca_c*D1/B2=0
B3=((2*alpha+k1_2ca_o+k2_2ca_o).*(1-A2)+3*beta+(k1_ca_o+k2_ca_o).*(1-A1).*A2+k2_ca_c.*A1.*A2);
D3=(nchoosek(4,1)*m_inf.*(1-m_inf).^3.*k1_ca_o+nchoosek(4,2)*m_inf.^2.*(1-m_inf).^2.*k1_2ca_o-k2_ca_c.*A1.*D2-(k1_ca_o+k2_ca_o).*(1-A1).*D2+(2*alpha+k1_2ca_o+k2_2ca_o).*D2)./B3;
A3=3*beta./B3;
B4=((alpha+k1_3ca_o+k2_3ca_o).*(1-A3)+4*beta+(k1_2ca_o+k2_2ca_o).*(1-A2).*A3+(k1_ca_o+k2_ca_o).*(1-A1).*A2.*A3+k2_ca_c.*A1.*A2.*A3);
D4=(nchoosek(4,1)*m_inf.*(1-m_inf).^3.*k1_ca_o+nchoosek(4,2)*m_inf.^2.*(1-m_inf).^2.*k1_2ca_o+nchoosek(4,3)*m_inf.^3.*(1-m_inf).*k1_3ca_o-k2_ca_c.*A1.*(D2+A2.*D3)-(k1_ca_o+k2_ca_o).*(1-A1).*(D2+A2.*D3)-(k1_2ca_o+k2_2ca_o).*((1-A2).*D3-D2)+(alpha+k1_3ca_o+k2_3ca_o).*D3)./B4;
A4=4*beta./B4;

%   ccccy=A1*A2*A3*A4*y+A1*A2*A3*D4+A1*A2*D3+A1*D2+D1;
%     cccoy=(1-A1)*A2*A3*A4*y+D2*(1-A1)+A2*D3*(1-A1)+(1-A1)*A2*A3*D4-D1;
%     ccooy=(1-A2)*A3*A4*y+A3*D4*(1-A2)+D3*(1-A2)-D2;
%     coooy=(1-A3)*A4*y+D4*(1-A3)-D3;
%     ooooy=(1-A4)*y-D4;

ca4_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
ca4_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
ca4_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);
ca4_o3_v=nchoosek(n_ca,3).*m_inf_v.^3.*(1-m_inf_v).^(n_ca-3);
ca4_o4_v=nchoosek(n_ca,4).*m_inf_v.^4.*(1-m_inf_v).^(n_ca-4);

tauBK_qss_4=1./(k2_ca_c.*A1.*A2.*A3.*A4+(k1_ca_o+k2_ca_o).*(1-A1).*A2.*A3.*A4 +(k1_2ca_o+k2_2ca_o).*(1-A2).*A3.*A4+(k1_3ca_o+k2_3ca_o).*(1-A3).*A4+(k1_4ca_o+k2_4ca_o).*(1-A4));
mBKinf_qss_4=(ca4_o1_v.*k1_ca_o+k1_2ca_o.*ca4_o2_v+k1_3ca_o.*ca4_o3_v+k1_4ca_o.*ca4_o4_v - k2_ca_c.*(D2.*A1+ D3.*A1.*A2+D4.*A1.*A2.*A3) ...
    -(k1_ca_o+k2_ca_o).*(1-A1).*(D2+A2.*D3+A2.*A3.*D4)-(k1_2ca_o+k2_2ca_o).*((1-A2).*(A3.*D4+D3)-D2)-(k1_3ca_o+k2_3ca_o).*((1-A3).*D4-D3)+(k1_4ca_o+k2_4ca_o).*D4).*tauBK_qss_4;

% tau (tauBK_qss_minf_4) and ss (mBKinf_qss_minf_4) BK activation given by  
% equation (29) with k=n_ca=4;

tauBK_qss_minf_4=1./((1-m_inf).^n_ca.*k2_ca_c+ca4_o1_v.*(k1_ca_o+k2_ca_o)+ca4_o2_v.*(k1_2ca_o+k2_2ca_o)+ca4_o3_v.*(k1_3ca_o+k2_3ca_o)+ca4_o4_v.*(k1_4ca_o+k2_4ca_o));
mBKinf_qss_minf_4=(ca4_o1_v.*k1_ca_o+ca4_o2_v.*k1_2ca_o+ca4_o3_v.*k1_3ca_o+ca4_o4_v.*k1_4ca_o).*tauBK_qss_minf_4;


f=figure(1);
hold on, 
grid on,
plot(V_vect,mBKinf_qss_1,'-c');
plot(V_vect,mBKinf_qss_minf_1,'--c');
plot(V_vect,mBKinf_qss_2,'-g');
plot(V_vect,mBKinf_qss_minf_2,'--g');
plot(V_vect,mBKinf_qss_4,'-r');
plot(V_vect,mBKinf_qss_minf_4,'--r');
xlim([-58 58])
plot(V_vect,m_inf,'--','color',[0.5 0.5 0.5]);
% set(gca,'xTickLabel',[])
% set(f,'position',[20 20 150 90])
% set(findall(f,'type','text'),'fontSize',8)
set(f,'position',[20 20 300 300])
set(gca,'fontSize',12);
title('Figure 2B (upper)')
xlabel(' V [mV])')
ylabel('m_{BK}, m_{CaV}') 
f=figure(2); 
hold on
grid on

plot(V_vect,tauBK_qss_1,'-c');
plot(V_vect,tauBK_qss_minf_1,'--c');
plot(V_vect,tauBK_qss_2,'-g');
plot(V_vect,tauBK_qss_minf_2,'--g');
plot(V_vect,tauBK_qss_4,'-r');
plot(V_vect,tauBK_qss_minf_4,'--r');
% plot(V_vect,tau_v,'--','color',[0.5 0.5 0.5]);
xlim([-58 58])
ylim([0 2])
% set(gca,'fontSize',8);
% set(f,'position',[20 20 150 60])
set(gca,'fontSize',12);
set(f,'position',[20 20 300 300])
title('Figure 2B (lower)')
xlabel(' V [mV])') 
ylabel(' \tau_{BK} [ms]') 