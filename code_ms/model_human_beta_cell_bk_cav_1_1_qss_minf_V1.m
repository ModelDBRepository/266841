

function xdot=model_human_beta_cell_bk_cav_1_1_qss_minf(t,x,flag)

% Human beta-cell model with 1:1 BK-CaV complexes. The I_BK current is modeled by
% equation (27) (with inactivating L- and T-type CaVs), or by equation (28) (with 
% non-inactivating P/Q-type CaVs), where the BK activation is modeled by equation (29)
% (assuming instantaneous activation of CaVs).

global global_cBK_T global_cBK_L global_cBK_PQ  nmCaPQ
global gBK gCaPQ VmCaPQ

global gkatp gleak gCaL gNa r_bk

%state variables

v=x(1);
mkv=x(2);
hNa=x(3);
hCaL=x(4);
hCaT=x(5);
xERG=x(6);
yERG=x(7);
mBK_L=x(8);
mBK_PQ=x(9);
mBK_T=x(10);

%%%%% current parameters %%%%%%%

%leak
%gleak=0.015;
vleak=-30;

%IKv
taumkv0=2;
Vmkv=0;
nmkv=-10;
gkv=1;

%IBK
%taumBK=2;
% VmBK=0;
% nmBK=-10;
% BBK=20;

% %IKir
% VmKir=-90;
% nmKir=15;
% gKir=0.1;

%Na current
% gNa=0.4;
VmNa=-18;
nmNa=-5;
VhNa=-42;
nhNa=6;
tauhNa=2;

%L-type Ca2+ current
%gCaL=0.14;
VmCaL=-25;
nmCaL=-6;
tauhCaL=20;

%PQ-type Ca2+ current
% gCaPQ=0.17;
% VmCaPQ=-10;
% nmCaPQ=-10;

%T-type Ca2+ current
gCaT=0.05;
VmCaT=-40;
nmCaT=-4;
VhCaT=-64;
nhCaT=8;
tauhCaT=7;

% %ISK current
% gSK=0.1;
% kSK=0.57;
% nSK=5.2;

%Katp current
% gkatp=0.015;
%gkatp=0.015;

%Herg current
VxERG=-30; 
nxERG=-10; 
tauxERG=100;
 
tauyERG=50; 
VyERG=-42; 
nyERG=17.5; 
gERG=0.2;


%Nernst voltages
VNa=70;
VCa=65;
VK=-75;
VCl=-40;


%%% Leak current
Ileak=gleak*(v-vleak);


%%% IKv
mkvinf=1/(1+exp((v-Vmkv)/nmkv));
taumkv=taumkv0+10*exp(min(log(3),(-20-v)/6));
IKv=gkv*mkv*(v-VK);
 

% %%% IKir
% IKir=gKir*(v-VK)/(1+exp((v-VmKir)/nmKir));


%%% Na current
hNainf=1/( 1+exp((v-VhNa)/nhNa) );
mNainf=1/( 1+exp((v-VmNa)/nmNa) );
INa=gNa*mNainf*hNa*(v-VNa);


%%% L-type Ca current
mCaLinf=1/( 1+exp((v-VmCaL)/nmCaL) );
hCaLinf=max(0,min(1,1+mCaLinf*(v-VCa)/57));
ICaL = gCaL*mCaLinf*hCaL*(v-VCa);



%%% PQ-type Ca current
mCaPQinf=1/( 1+exp((v-VmCaPQ)/nmCaPQ) );
ICaPQ = gCaPQ*mCaPQinf*(v-VCa);


%%% T-type Ca current
mCaTinf=1/( 1+exp((v-VmCaT)/nmCaT) );
hCaTinf=1/( 1+exp((v-VhCaT)/nhCaT) );
ICaT = gCaT*mCaTinf*hCaT*(v-VCa);


%%% hERG-channels
%%% Based on Rosati et al. 2000 , tau-act from Schonherr 1999


xERGinf=1/( 1+exp((v-VxERG)/nxERG) );
yERGinf=1/( 1+exp((v-VyERG)/nyERG) );

IERG = gERG*xERG*yERG*(v-VK);


%%% Parameters for the IBK current
    
%%% compute calcium level at 7nm of the pore (cav sensor) and ....
D_ca=250; %microm^2 s^-1
F=9.6485*10^4; %C mol^-1
conv_F=10^(-15);  % mol a M/microm^3
conv_microM=10^6; % M to microM
%%%%r_ca=7*10^(-3); %microm
k_B=500; %microM^-1 * s^-1
B_tot=30; %microM
Eca=65; %mV
conv_V=10^(-3); %mV to V
g_ca_PQ=2.8; %pS 
g_ca_L=2.2; %pS 
g_ca_T=1.7; %pS 

conv_S=10^(-12); %pS to S
% Ca_rca_lev=i_ca/(8*pi*D_ca*F*conv_F*r_ca)*exp(-r_ca/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
%%%...and at the bk sensor r_bk
%%%%r_bk=30*10^(-3);% microm

i_ca_L=abs((v-Eca))*conv_V*g_ca_L*conv_S; % C/sec
i_ca_PQ=abs((v-Eca))*conv_V*g_ca_PQ*conv_S; % C/sec
i_ca_T=abs((v-Eca))*conv_V*g_ca_T*conv_S; % C/sec

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

k1=k1_0*exp(-alpha0*v);
k2=k2_0*exp(-beta0*v);


ca_c=0.2;
% ca_o=n_ca_o*(-6/20*v+ca_o_v0);
fco_Ca_c=k1*(ca_c^n1)/(ca_c^n1+K1^n1);
foc_Ca_c=k2*(K2^n2)/(K2^n2+ca_c^n2);

k1_ca_c=fco_Ca_c;
k2_ca_c=foc_Ca_c;

%%% IBK
IBK=gBK*(global_cBK_L*hCaL*mBK_L+global_cBK_PQ*mBK_PQ+global_cBK_T*hCaT*mBK_T)*(v-VK);

%%% Katp current
IKatp=gkatp*(v-VK);


%ODEs for BK activation

%%%%%%%%%%%%%%%%%%%
% BK activated by L-type %
%%%%%%%%%%%%%%%%%%%

Ca_rbk_lev_L=i_ca_L/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_o=Ca_rbk_lev_L;

fco_Ca_o=k1*(ca_o^n1)/(ca_o^n1+K1^n1);
foc_Ca_o=k2*(K2^n2)/(K2^n2+ca_o^n2);
k1_ca_o=fco_Ca_o;
k2_ca_o=foc_Ca_o;

y=mBK_L;
m_inf=mCaLinf;
cy=(1-m_inf)*y;       
oy=m_inf*y;       
mBK_L_dot=-k2_ca_c*cy+k1_ca_o*(m_inf-oy)-k2_ca_o*oy;

%%%%%%%%%%%%%%%%%%%%%
% BK activated by P/Q-type %
%%%%%%%%%%%%%%%%%%%%%

Ca_rbk_lev_PQ=i_ca_PQ/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_o=Ca_rbk_lev_PQ;

fco_Ca_o=k1*(ca_o^n1)/(ca_o^n1+K1^n1);
foc_Ca_o=k2*(K2^n2)/(K2^n2+ca_o^n2);
k1_ca_o=fco_Ca_o;
k2_ca_o=foc_Ca_o;

y=mBK_PQ;
m_inf=mCaPQinf;
cy=(1-m_inf)*y;       
oy=m_inf*y;       
mBK_PQ_dot=-k2_ca_c*cy+k1_ca_o*(m_inf-oy)-k2_ca_o*oy;


%%%%%%%%%%%%%%%%%%%
% BK activated by T-type %
%%%%%%%%%%%%%%%%%%%

Ca_rbk_lev_T=i_ca_T/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_o=Ca_rbk_lev_T;

fco_Ca_o=k1*(ca_o^n1)/(ca_o^n1+K1^n1);
foc_Ca_o=k2*(K2^n2)/(K2^n2+ca_o^n2);
k1_ca_o=fco_Ca_o;
k2_ca_o=foc_Ca_o;

y=mBK_T;
m_inf=mCaTinf;
cy=(1-m_inf)*y;       
oy=m_inf*y;       
mBK_T_dot=-k2_ca_c*cy+k1_ca_o*(m_inf-oy)-k2_ca_o*oy;

% other ODEs

v_dot= -(IERG + IBK + IKv + INa + ICaL + ICaPQ + ICaT + IKatp + Ileak );

mkv_dot=(mkvinf-mkv)/taumkv;

% mBK_dot=(mBKinf-mBK)/taumBK;

hNa_dot=(hNainf-hNa)/tauhNa;

hCaL_dot=(hCaLinf-hCaL)/tauhCaL;

hCaT_dot=(hCaTinf-hCaT)/tauhCaT;


%activation ERG
xERG_dot= (xERGinf-xERG)/tauxERG;
%inactivation ERG
yERG_dot = (yERGinf-yERG)/tauyERG;

xdot=[v_dot;mkv_dot;hNa_dot;hCaL_dot;hCaT_dot;xERG_dot;yERG_dot;mBK_L_dot;mBK_PQ_dot;mBK_T_dot];