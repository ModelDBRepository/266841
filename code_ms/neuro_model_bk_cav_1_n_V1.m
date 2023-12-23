function x_dot=neuro_model_bk_cav_1_n(t,x,flag,tu,u_in)

% Neuronal model with 1:n BK-CaV complexes. The I_BK current is modeled by
% equation (28) (BK coupled with non-inactivating  CaVs), 
% where the BK activation is modeled by equation (26)
% (assuming CaV kinetics), or by equation (29) (assuming instantaneous
% activation of CaVs, tag_minf=1 in the code)

global gBK no_leak_na_k gKA gKAHP gKSOR n_ca_o global_cBK_L r_bk  VmCa n_ind tau_inf tag_minf kmCa

%state variables
v=x(1);
hNa=x(2);
mKDR=x(3);
mKA=x(4);
hKA=x(5);
mKSOR=x(6);
mBK_L=x(7);
Ca_BK=x(8);
Ca_SK=x(9);
Ca_i=x(10);
a=x(11);
mCa=x(12);

u_0=0;
if isempty(tu)
    u=u_0;
else
    
    if length(tu)==1
        if t<tu
            u=u_0;       
        else
            u=u_in;
        end    
    else
        if t<tu(1)
            u=u_0;
        elseif t>=tu(end)
            u=u_in(end);
        else
            k=1;
            t_found=false;
            while not(t_found)
                if t>=tu(k) && t<tu(k+1)
                    t_found=true;
                    u=u_in(k);
                end
                k=k+1;
            end
        end
    end
end


Ca_r=113; %nM
Ca_o=4e6; %nM


% Reverse potentials
E_Na=50;
E_K=-96; 
E_leak=-65;
E_Ca=12.5*log(Ca_o/Ca_i);

%current parameters

% Na
gNa=14;
VmNa=38;
kmNa=4;
VhNa=45;
khNa=2;

minfNa=1/(1+exp((-v-VmNa)/kmNa));
hinfNa=1/(1+exp((v+VhNa)/khNa));

tauhNa=20/((1+exp((v-1)/14))*(1+exp((v+30)/34)))- 0.075;

% Ca
gCa=0.1;
% VmCa=10;
%kmCa=7;

minfCa=1/(1+exp((-v-VmCa)/kmCa));

% K-DR
gKDR=14;
VmKDR=-2;
kmKDR=11;

minfKDR=1/(1+exp((-v-VmKDR)/kmKDR));

taumKDR=8/(exp((v+40)/19)+exp((-v-40)/20))+1.8;

% K-A
%gKA=14;
VmKA=45;
kmKA=11;
VhKA=80;
khKA=6.5;

minfKA=1/(1+exp((-v-VmKA)/kmKA));
hinfKA=1/(1+exp((v+VhKA)/khKA));

taumKA=35/(1.1*exp((v+43)/17)+exp((-v-20)/13));
tauhKA=5+120/(exp((v+75)/30)+exp((-v-45)/18));

%K-SOR 
% gKSOR=0.06;
VmKSOR=60;
kmKSOR=1.8;

minfKSOR=1/(1+exp((-v-VmKSOR)/kmKSOR));
taumKSOR=260;

%K-AHP
% gKAHP=0.18;
qinf=1/(1+exp(-1.12-2.508*log((Ca_SK-Ca_r)/1000)));

%K-BK
%pinf=1/((1+470/(Ca_BK^2.38))*(1+exp((-v-140*log10(Ca_BK)+370)/7.4)));

%%% compute calcium level at 7nm of the pore (cav sensor) and at the bk
%%% sensor (r_bk)
D_ca=250; %microm^2 s^-1
F=9.6485*10^4; %C mol^-1
conv_F=10^(-15);  % mol a M/microm^3
conv_microM=10^6; % M to microM
%%%%r_ca=7*10^(-3); %microm
k_B=500; %microM^-1 * s^-1
B_tot=30; %microM
% Eca=65; %mV
conv_V=10^(-3); %mV to V
g_ca_PQ=2.8; %pS 
g_ca_L=2.2; %pS 
g_ca_T=1.7; %pS 
conv_S=10^(-12); %pS to S
% Ca_rca_lev=i_ca/(8*pi*D_ca*F*conv_F*r_ca)*exp(-r_ca/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
%%%...and at the bk sensor r_bk
%%%%r_bk=30*10^(-3);% microm

i_ca_L=abs((v-E_Ca))*conv_V*g_ca_L*conv_S; % C/sec
% i_ca_PQ=abs((v-E_Ca))*conv_V*g_ca_PQ*conv_S; % C/sec
% i_ca_T=abs((v-E_Ca))*conv_V*g_ca_T*conv_S; % C/sec

%%%% BK parameters
BK_par=[  1.1093   3.3206  2.3298    0.0223    1.6012   16.5793    0.1000    0.4614];

k1_0=BK_par(1);
k2_0=BK_par(2);
K1=BK_par(6);
K2=BK_par(7);
n1=BK_par(3);
n2=BK_par(8);
beta0=BK_par(4);
alpha0=-BK_par(5)*beta0;
k1=k1_0*exp(-alpha0*v);
k2=k2_0*exp(-beta0*v);

ca_c=0.2;
% ca_o=n_ca_o*(-6/20*v+ca_o_v0);
fco_Ca_c=k1*(ca_c^n1)/(ca_c^n1+K1^n1);
foc_Ca_c=k2*(K2^n2)/(K2^n2+ca_c^n2);

k1_ca_c=fco_Ca_c;
k2_ca_c=foc_Ca_c;

%%% L-type Ca2+ channel
mCaLinf=minfCa;
betaL= (1-mCaLinf)/tau_inf;    
alphaL=mCaLinf/tau_inf;

%%%  n_ca_o is the number of  synchronized CaV channels (assumed n_ca_o=1)
alpha=n_ca_o*alphaL;
beta=n_ca_o*betaL;

%%% calcium level at the bk sensor (r_bk)
Ca_rbk_lev_L=i_ca_L/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_o=n_ca_o*Ca_rbk_lev_L;

fco_Ca_o=k1*(ca_o^n1)/(ca_o^n1+K1^n1);
foc_Ca_o=k2*(K2^n2)/(K2^n2+ca_o^n2);
k1_ca_o=fco_Ca_o;
k2_ca_o=foc_Ca_o;

fco_2Ca_o=k1*((2*ca_o)^n1)/((2*ca_o)^n1+K1^n1);
foc_2Ca_o=k2*(K2^n2)/(K2^n2+(2*ca_o)^n2);

k1_2ca_o=fco_2Ca_o;
k2_2ca_o=foc_2Ca_o;

fco_3Ca_o=k1*((3*ca_o)^n1)/((3*ca_o)^n1+K1^n1);
foc_3Ca_o=k2*(K2^n2)/(K2^n2+(3*ca_o)^n2);

k1_3ca_o=fco_3Ca_o;
k2_3ca_o=foc_3Ca_o;

fco_4Ca_o=k1*((4*ca_o)^n1)/((4*ca_o)^n1+K1^n1);
foc_4Ca_o=k2*(K2^n2)/(K2^n2+(4*ca_o)^n2);

k1_4ca_o=fco_4Ca_o;
k2_4ca_o=foc_4Ca_o;

y=mBK_L; % state variable for BK activation
if n_ind==1 % number of independent CaV channels
    
    if tag_minf % if equal to 1 then assume that oy=m_inf*y 
            
        m_inf=minfCa;
        cy=(1-m_inf)*y;
        oy=m_inf*y;
        mBK_L_dot=-k2_ca_c*cy+k1_ca_o*(mCa-oy)-k2_ca_o*oy; % ode for BK activation
    else
        
        m_inf=mCa;
        tau_BK=(alpha+beta+k2_ca_c)/((alpha+k2_ca_c)*(k1_ca_o+k2_ca_o)+beta*k2_ca_c);
        mBK_inf=k1_ca_o*m_inf*tau_BK;
        
        mBK_L_dot=(mBK_inf-y)/tau_BK; % ode for BK activation
    end


elseif n_ind==2
    
    if tag_minf % instantaneous activation of CaVs
        n_ca=2; 
        m_inf=minfCa;   
        m_inf_v=m_inf;
        ca2_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
        ca2_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
        ca2_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);
     
        ccy=ca2_o0_v*y;
        coy=ca2_o1_v*y;
        ooy=ca2_o2_v*y;
    else
        m_inf=mCa;
        A=beta/(2*alpha+k2_ca_c+beta);  % A_1 in the paper
        B_0=((alpha+k1_ca_o+k2_ca_o)*(1-A)+2*beta+k2_ca_c*A);  % B_2 in the paper
        B=2*beta/B_0;  % A_2 in the paper
        D=2*(1-m_inf)*m_inf*k1_ca_o/B_0;  %D_2 in the paper


        ooy=(1-B)*y-D;
        coy=(1-A)*(B*y+D);
        ccy=A*(B*y+D);
    end
    
    % ode for BK activation
    mBK_L_dot=-k2_ca_c*ccy+k1_ca_o*(2*mCa*(1-mCa)-coy)-k2_ca_o*coy-k2_2ca_o*ooy+k1_2ca_o*(mCa^2-ooy);
    
elseif n_ind==3 % instantaneous activation of CaVs
    
    if tag_minf
        
        n_ca=3;
        m_inf=minfCa;
        m_inf_v=m_inf;
        ca3_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
        ca3_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
        ca3_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);
        ca3_o3_v=nchoosek(n_ca,3).*m_inf_v.^3.*(1-m_inf_v).^(n_ca-3);
     
        cccy=ca3_o0_v*y;
        ccoy=ca3_o1_v*y;
        cooy=ca3_o2_v*y;
        oooy=ca3_o3_v*y;
        
    else
        m_inf=mCa;
        
        A1=beta/(3*alpha+k2_ca_c+beta);
        D1=0;
        B2=((2*alpha+k1_ca_o+k2_ca_o)*(1-A1)+2*beta+k2_ca_c*A1);
        A2=2*beta/B2;
        D2=3*m_inf*(1-m_inf)^2*k1_ca_o/B2; %-k2_ca_c*D1/B2=0
        B3=((alpha+k1_2ca_o+k2_2ca_o)*(1-A2)+3*beta+(k1_ca_o+k2_ca_o)*(1-A1)*A2+k2_ca_c*A1*A2);
        D3=(nchoosek(3,1)*m_inf*(1-m_inf)^2*k1_ca_o+nchoosek(3,2)*m_inf^2*(1-m_inf)*k1_2ca_o-k2_ca_c*A1*D2-(k1_ca_o+k2_ca_o)*(1-A1)*D2+(alpha+k1_2ca_o+k2_2ca_o)*D2)/B3;
        A3=3*beta/B3;

        cccy=A1*A2*A3*y+A1*A2*D3+A1*D2+D1;
        ccoy=(1-A1)*A2*A3*y+D2*(1-A1)+A2*D3*(1-A1)-D1;
        cooy=(1-A2)*A3*y+D3*(1-A2)-D2;
        oooy=(1-A3)*y-D3;
    end
    
    % ode for BK activation
    mBK_L_dot=-k2_ca_c*cccy+k1_ca_o*(3*m_inf*(1-m_inf)^2-ccoy)-k2_ca_o*ccoy-k2_2ca_o*cooy+k1_2ca_o*(3*m_inf^2*(1-m_inf)-cooy)-k2_3ca_o*oooy+k1_3ca_o*(m_inf^3-oooy);
   
elseif n_ind==4 % instantaneous activation of CaVs
     
    if tag_minf
        n_ca=4;
        m_inf=minfCa;
        m_inf_v=m_inf;
        ca4_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
        ca4_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
        ca4_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);
        ca4_o3_v=nchoosek(n_ca,3).*m_inf_v.^3.*(1-m_inf_v).^(n_ca-3);
        ca4_o4_v=nchoosek(n_ca,4).*m_inf_v.^4.*(1-m_inf_v).^(n_ca-4);

        
        ccccy=ca4_o0_v*y;
        cccoy=ca4_o1_v*y;
        ccooy=ca4_o2_v*y;
        coooy=ca4_o3_v*y;
        ooooy=ca4_o4_v*y;
        
    else
        
        m_inf=mCa;
        A1=beta/(4*alpha+k2_ca_c+beta);
        D1=0;
        B2=((3*alpha+k1_ca_o+k2_ca_o)*(1-A1)+2*beta+k2_ca_c*A1);
        A2=2*beta/B2;
        D2=4*m_inf*(1-m_inf)^3*k1_ca_o/B2; % k2_ca_c*D1/B2=0
        B3=((2*alpha+k1_2ca_o+k2_2ca_o)*(1-A2)+3*beta+(k1_ca_o+k2_ca_o)*(1-A1)*A2+k2_ca_c*A1*A2);
        D3=(nchoosek(4,1)*m_inf*(1-m_inf)^3*k1_ca_o+nchoosek(4,2)*m_inf^2*(1-m_inf)^2*k1_2ca_o-k2_ca_c*A1*D2-(k1_ca_o+k2_ca_o)*(1-A1)*D2+(2*alpha+k1_2ca_o+k2_2ca_o)*D2)/B3;
        A3=3*beta/B3;
        B4=((alpha+k1_3ca_o+k2_3ca_o)*(1-A3)+4*beta+(k1_2ca_o+k2_2ca_o)*(1-A2)*A3+(k1_ca_o+k2_ca_o)*(1-A1)*A2*A3+k2_ca_c*A1*A2*A3);
        D4=(nchoosek(4,1)*m_inf*(1-m_inf)^3*k1_ca_o+nchoosek(4,2)*m_inf^2*(1-m_inf)^2*k1_2ca_o+nchoosek(4,3)*m_inf^3*(1-m_inf)*k1_3ca_o-k2_ca_c*A1*(D2+A2*D3)-(k1_ca_o+k2_ca_o)*(1-A1)*(D2+A2*D3)-(k1_2ca_o+k2_2ca_o)*((1-A2)*D3-D2)+(alpha+k1_3ca_o+k2_3ca_o)*D3)/B4;
        A4=4*beta/B4;     
        
        ccccy=A1*A2*A3*A4*y+A1*A2*A3*D4+A1*A2*D3+A1*D2+D1;
        cccoy=(1-A1)*A2*A3*A4*y+D2*(1-A1)+A2*D3*(1-A1)+(1-A1)*A2*A3*D4-D1;
        ccooy=(1-A2)*A3*A4*y+A3*D4*(1-A2)+D3*(1-A2)-D2;
        coooy=(1-A3)*A4*y+D4*(1-A3)-D3;
        ooooy=(1-A4)*y-D4;
                
    end
    % ode for BK activation
    mBK_L_dot=-k2_ca_c*ccccy+k1_ca_o*(4*m_inf*(1-m_inf)^3-cccoy)-k2_ca_o*cccoy-k2_2ca_o*ccooy+k1_2ca_o*(6*m_inf^2*(1-m_inf)^2-ccooy)-k2_3ca_o*coooy+k1_3ca_o*(4*m_inf^3*(1-m_inf)-coooy)-k2_4ca_o*ooooy+k1_4ca_o*(m_inf^4-ooooy);

   
end

%%% I_BK current
iBK=gBK*(global_cBK_L*mBK_L)*(v-E_K);
%iBK=gBK*p*(v-E_K);

% leak
gleak=0.083;
gNa_leak=0.018; %
gK_leak=0.066;%

ka=50;
kb=8;
Vb=-60;
tau_a=75;
ainf=tanh((Ca_i-Ca_r)/ka);
binf=1/(1+exp(-(v+Vb)/kb));

% other currents
iNa=gNa*minfNa^3*hNa*(v-E_Na);
iCa=gCa*minfCa^1*(v-E_Ca);
iKDR=gKDR*mKDR^3*(v-E_K);
iKA=gKA*mKA^4*hKA*(v-E_K);
iKSOR=gKSOR*mKSOR*(v-E_K);
iKAHP=gKAHP*qinf^2*(v-E_K);

if no_leak_na_k
    i_leak=gleak*(v-E_leak);
else
    iNa_leak=gNa_leak*(v-E_Na);

    f_leak=a*binf;
    iK_leak=gK_leak*(1-f_leak)*(v-E_K);
    i_leak=iK_leak+iNa_leak;
end

%ODEs

alpha_BK=100;
alpha_SK=1.6;
alpha_i=1.4;
tau_BK=1;
tau_SK=656;
tau_i=2.33*1e3;
C=1;

v_dot=-1/C*(iNa+iCa+iKDR+iKA+iBK+iKAHP+iKSOR+i_leak+u);

Ca_i_dot=-alpha_i*iCa-1/tau_i*(Ca_i-Ca_r);

hNa_dot=(hinfNa-hNa)/tauhNa;

mKDR_dot=(minfKDR-mKDR)/taumKDR;

mKA_dot=(minfKA-mKA)/taumKA;

hKA_dot=(hinfKA-hKA)/tauhKA;

mKSOR_dot=(minfKSOR-mKSOR)/taumKSOR;

% p_dot=(pinf-p)/tau_p;

CaBK_dot=-alpha_BK*iCa-1/tau_BK*(Ca_BK-Ca_r);

CaSK_dot=-alpha_SK*iCa-1/tau_SK*(Ca_SK-Ca_r);

a_dot=(ainf-a)/tau_a;

mCa_dot=(minfCa-mCa)*n_ca_o/tau_inf;

x_dot=[v_dot;hNa_dot;mKDR_dot;mKA_dot;hKA_dot;mKSOR_dot;mBK_L_dot;CaBK_dot;CaSK_dot;Ca_i_dot;a_dot;mCa_dot];
