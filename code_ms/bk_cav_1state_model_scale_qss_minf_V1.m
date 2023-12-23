
function xdot=bk_cav_1state_model_scale_qss_minf(t,x,flag,tv,V_in)

% 1-ODE model for the BK activation (Eq. 29 with k=1) with m_CaV= m_{CaV,\infty}  in 1:1
% BK-CaV complex

global k1_0 K1  k2_0 K2 n1 n2 alpha0 beta0 V_0 ca_c beta1  alpha1 beta1_0 alpha1_0 rho

y=x(1);
%c=x(2);

if isempty(tv)
    V=V_0;
else
    
    if length(tv)==1
        if t<tv
            V=V_0;       
        else
            V=V_in;
        end    
    else
        if t<tv(1)
            V=V_0;
        elseif t>tv(end)
            V=V_in(end);
        else
            k=1;
            t_found=false;
            while not(t_found)
                if t>=tv(k) && t<=tv(k+1)
                    t_found=true;
                    V=V_in(k);
                end
                k=k+1;
            end
        end
    end
end

% BK-CaV parameters

alpha=alpha1_0*exp(-alpha1*V); 
beta_0=beta1_0*exp(-beta1*V);
beta=(beta_0+alpha)* rho;

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


k1=k1_0*exp(-alpha0*V);

k2=k2_0*exp(-beta0*V);

%fco_Ca_c=(ca_c^n1)/(ca_c^n1+K1^n1);
fco_Ca_o=(ca_o^n1)/(ca_o^n1+K1^n1);

foc_Ca_c=(K2^n2)/(K2^n2+ca_c^n2);
foc_Ca_o=(K2^n2)/(K2^n2+ca_o^n2);

k1_ca_o=k1*fco_Ca_o;
%k1_ca_c=k1*fco_Ca_c;
k2_ca_o=k2*foc_Ca_o;
k2_ca_c=k2*foc_Ca_c;

% o=1-c;
% dot_c=-alpha*c+beta*o;

m_inf=alpha/(alpha+beta);

%m_inf=o;


m_inf_v=m_inf;
n_ca=1;
ca1_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
ca1_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);

cy=ca1_o0_v*y;
oy=ca1_o1_v*y;


%1-ODE
ydot=-k2_ca_c*cy+k1_ca_o*(m_inf-oy)-k2_ca_o*oy;


xdot=ydot;