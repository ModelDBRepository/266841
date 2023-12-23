
function zdot=bk_n_cav_states_model_scale(t,z,flag,tv,V_in)

% 2*(n+1)-state ODE model for the BK activation in 1:n BK-CaV complex with
% n non-inactivated CaVs (n=n_c) (Eqs. S19-S24)

global k1_0 K1  k2_0 K2 n1 n2 alpha0 beta0 V_0 ca_c n_c beta1  alpha1 beta1_0 alpha1_0 rho

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

xp_dot=zeros(1,n_c);
yp_dot=zeros(1,n_c);

for p=0:n_c
    px=p+1;
    py=p+1+(n_c+1);
    
    if p==0
        p0x=0;
        p0y=0;
        ca_b=ca_c;
    else    
        p0x=(n_c-(p-1))*alpha*z(px-1);
        p0y=(n_c-(p-1))*alpha*z(py-1);
        ca_b=ca_o*p;
    end
    
    fco_Ca=(ca_b^n1)/(ca_b^n1+K1^n1);
    foc_Ca=(K2^n2)/(K2^n2+ca_b^n2);

    
    if p==n_c
        pnx=0;
        pny=0;
    else
        pnx=(p+1)*beta*z(px+1);
        pny=(p+1)*beta*z(py+1);
    end
    
    % C_{n-p}O_p where n is the number of calcium channels and p the number
    % of open channels: for example n=2 p=0, C_2_O_0=>CC, n=2,p=1=>CO
    
    % C_{n-p}O_pXdot=-(n-p)*alpha*C_{n-p}O_pX
    % +(n-(p-1))*alpha*C_{n-(p-1)}*O_{p-1}X+
    % (p+1)*beta*C_{n-(p+1)}*O_{p+1}X -p*beta*C_{n-p}*O_pX +
    % k2*foc(Ca(O_p))*C_{n-p}O_pY -k1*fco(Ca(O_p))*C_{n-p}O_pX
    
    xp_dot(px)=-(n_c-p)*alpha*z(px)+ p0x+pnx-p*beta*z(px) + k2*foc_Ca*z(py) -k1*fco_Ca*z(px);
     
    yp_dot(px)=-(n_c-p)*alpha*z(py)+ p0y+pny-p*beta*z(py) - k2*foc_Ca*z(py) + k1*fco_Ca*z(px);

end

zdot=[xp_dot'; yp_dot'];