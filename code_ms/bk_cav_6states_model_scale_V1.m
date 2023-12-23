
function xdot=bk_cav_6states_model_scale(t,x,flag,tv,V_in)

% 6-state ODE model for 1:1 BK-CaV with inactivating CaV (Eqs. S6-S11)

global k1_0 K1  k2_0 K2 n1 n2 alpha0 beta0 V_0  ca_c N_tot k_minus k_plus beta1 alpha1 beta1_0 alpha1_0 rho

oy=x(1);

cy=x(2);

ox=x(3);

bx=x(4);

by=x(5);

cx=N_tot-oy-cy-ox-bx -by; 

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

alpha=alpha1_0*exp(-alpha1*V);
% beta_s=0.1650*exp(-0.1735*V);
beta_s=beta1_0*exp(-beta1*V);
beta=(beta_s+alpha)* rho;

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

fco_Ca_c=(ca_c^n1)/(ca_c^n1+K1^n1);
fco_Ca_o=(ca_o^n1)/(ca_o^n1+K1^n1);

foc_Ca_c=(K2^n2)/(K2^n2+ca_c^n2);
foc_Ca_o=(K2^n2)/(K2^n2+ca_o^n2);

dot_oy=k1*fco_Ca_o*ox -k2*foc_Ca_o*oy + alpha*cy - beta*oy + k_minus*by - k_plus*ca_o_rca*oy;

dot_cy=k1*fco_Ca_c*cx -k2*foc_Ca_c*cy - alpha*cy + beta*oy;

dot_ox= alpha*cx - beta*ox -k1*fco_Ca_o*ox + k2*foc_Ca_o*oy + k_minus*bx - k_plus*ca_o_rca*ox;

dot_bx= -k_minus*bx + k_plus*ca_o_rca*ox - k1*fco_Ca_c*bx + k2*foc_Ca_c*by;

dot_by= -k_minus*by + k_plus*ca_o_rca*oy + k1*fco_Ca_c*bx - k2*foc_Ca_c*by;

xdot=[dot_oy;dot_cy;dot_ox;dot_bx;dot_by];