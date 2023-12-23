function dy=model_pituitary_ODE(t,Y,flag,struct)
%% A pituitary lactotroph model with BK coupled with n CaVs. 
% The I_BK current is modeled by equation (28) (BK coupled with non-inactivating  CaVs), 
% where the BK activation is modeled by the complete ODE model with 2*(n+1) states 
% described by equations (S19)-(S24).
%% model parameters
C=struct.C;  %pF
gcal=struct.gca;  %nS
Vca=struct.Vca; %mV
Vm=struct.Vm; %mV
sm=struct.sm;  %mV
gk=struct.gk; %nS
Vk=struct.Vk; %mV
Vn=struct.vn;  %mV
sn=struct.sn;
taun=struct.taun;%ms
gl=struct.gl;  %nS
Vl=struct.Vl; %mV
ff=struct.fc;
alfa=struct.alfa; %microM*fC-1
kc=struct.Kc; %ms
gbk=struct.gbk;
gsk=struct.gsk;
ks=struct.ks;
chn=struct.nch;
num_ca=chn;
taum=1.25;

%% states
V=Y(1);
n=Y(2);  % open K 
c=Y(3);  % Ca



%% BK channels parameter

D_ca=250; %microm^2 s^-1
F=9.6485*10^4; %C mol^-1
conv_F=10^(-15);  % mol a M/microm^3
conv_microM=10^6; % M to microM
k_B=500; %microM^-1 * s^-1
B_tot=30; %microM
conv_V=10^(-3); %mV to V
g_ca_L=2; %pS 
conv_S=10^(-12); %pS to S
r_bk=30*10^(-3);% microm
i_ca_single=abs((V-Vca))*conv_V*g_ca_L*conv_S; % C/sec
bk_par=[  1.1093   3.3206  2.3298    0.0223    1.6012   16.5793    0.1000    0.4614];
k1_0=bk_par(1);
k2_0=bk_par(2);
K1=bk_par(6);
K2=bk_par(7);
n1=bk_par(3);
n2=bk_par(8);
beta0=bk_par(4);
alpha0=-bk_par(5)*beta0;
ca_c=0.2;
k1=k1_0*exp(-alpha0*V);
k2=k2_0*exp(-beta0*V);
foc_Ca_c=k2*(K2^n2)/(K2^n2+ca_c^n2);
k2_ca_c=foc_Ca_c;



%% Steady States
phik=1/(1+exp((Vn-V)/sn));
m_inf=1/(1+exp((Vm-V)/sm));
cinf=c^2/(c^2+ks^2);
alpha=m_inf/taum;
beta=(1-m_inf)/taum;

%% BK fluxes
ca_o=struct.nca*i_ca_single/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_c=0.2;
fco_Ca_o=k1*(ca_o^n1)/(ca_o^n1+K1^n1);
foc_Ca_o=k2*(K2^n2)/(K2^n2+ca_o^n2);
foc_Ca_c=k2*(K2^n2)/(K2^n2+ca_c^n2);
% k1_ca_o=fco_Ca_o;
% k2_ca_o=foc_Ca_o;
% k2_ca_c=foc_Ca_c;%k-C e

fco_2Ca_o=k1*((2*ca_o)^n1)/((2*ca_o)^n1+K1^n1);
foc_2Ca_o=k2*(K2^n2)/(K2^n2+(2*ca_o)^n2);

% k1_2ca_o=fco_2Ca_o;
% k2_2ca_o=foc_2Ca_o;

fco_3Ca_o=k1*((3*ca_o)^n1)/((3*ca_o)^n1+K1^n1);
foc_3Ca_o=k2*(K2^n2)/(K2^n2+(3*ca_o)^n2);

fco_4Ca_o=k1*((4*ca_o)^n1)/((4*ca_o)^n1+K1^n1);
foc_4Ca_o=k2*(K2^n2)/(K2^n2+(4*ca_o)^n2);

%% BK 1 Ca %%%%%%%%%
if chn==1

cy=1-Y(4)-Y(5)-Y(6);
oy=Y(6);
dY(4)= foc_Ca_c*cy+beta*Y(5)-alpha*Y(4);
dY(5)=alpha*Y(4)+foc_Ca_o*Y(6)-(beta+fco_Ca_o)*Y(5);
dY(6)=alpha*cy+fco_Ca_o*Y(5)-(beta+foc_Ca_o)*Y(6);
mBK=cy+oy;% open BK
end
%% BK 2 Ca %%chn%%%%%%%
if chn==2


ccy=1-Y(4)-Y(5)-Y(6)-Y(7)-Y(8);
coy=Y(8);
ooy=Y(7);
dY(4)=beta*Y(5)+foc_Ca_c*ccy-2*alpha*Y(4);
dY(5)=2*alpha*Y(4)+foc_Ca_o*Y(8)+2*beta*Y(6)-(beta+alpha+fco_Ca_o)*Y(5);
dY(6)=alpha*Y(5)+fco_2Ca_o*Y(7)-(2*beta+foc_2Ca_o)*Y(6);
dY(7)=foc_2Ca_o*Y(6)+alpha*Y(8)-(fco_2Ca_o+2*beta)*Y(7);
dY(8)=2*alpha*ccy+fco_Ca_o*Y(5)+2*beta*Y(7)-(beta+foc_Ca_o+alpha)*Y(8);
mBK=ccy+coy+ooy;
    
end
%% BK 3 Ca %%%%%%%%%
if chn==3

cccy=1-Y(4)-Y(5)-Y(6)-Y(7)-Y(8)-Y(9)-Y(10);
ccoy=Y(10);
cooy=Y(9);
oooy=Y(8);
dY(4)=beta*Y(5)+foc_Ca_c*cccy-3*alpha*Y(4);
dY(5)=3*alpha*Y(4)+foc_Ca_o*Y(10)+2*beta*Y(6)-(beta+2*alpha+fco_Ca_o)*Y(5);
dY(6)=2*alpha*Y(5)+foc_2Ca_o*Y(9)+3*beta*Y(7)-(2*beta+fco_2Ca_o+alpha)*Y(6);
dY(7)=alpha*Y(6)+foc_3Ca_o*Y(8)-(3*beta+fco_3Ca_o)*Y(7);
dY(8)=alpha*Y(9)+fco_3Ca_o*Y(7)-(3*beta+foc_3Ca_o)*Y(8);
dY(9)=fco_2Ca_o*Y(6)+2*alpha*Y(10)+3*beta*Y(8)-(foc_2Ca_o+2*beta+alpha)*Y(9);
dY(10)=3*alpha*cccy+fco_Ca_o*Y(5)+2*beta*Y(9)-(beta+foc_Ca_o+2*alpha)*Y(10);
mBK=cccy+ccoy+cooy+oooy;
end

%% BK 4 Ca %%%%%%%%%
if chn==4

ccccy =1-Y(4)-Y(5)-Y(6)-Y(7)-Y(8)-Y(9)-Y(10)-Y(11)-Y(12);
cccoy = Y(12);
ccooy = Y(11);
coooy = Y(10);
ooooy = Y(9);
dY(4)=beta*Y(5)+foc_Ca_c*ccccy-4*alpha*Y(4); 
dY(5)=4*alpha*Y(4)+foc_Ca_o*Y(12)+2*beta*Y(6)-(beta+3*alpha+fco_Ca_o)*Y(5);
dY(6)=3*alpha*Y(5)+foc_2Ca_o*Y(11)+3*beta*Y(7)-(2*beta+fco_2Ca_o+2*alpha)*Y(6);
dY(7)=2*alpha*Y(6)+foc_3Ca_o*Y(10)+4*beta*Y(8)-(3*beta+fco_3Ca_o+alpha)*Y(7);
dY(8)=alpha*Y(7)+foc_4Ca_o*Y(9)-(4*beta+fco_4Ca_o)*Y(8);
dY(9)=alpha*Y(10)+fco_4Ca_o*Y(8)-(4*beta+foc_4Ca_o)*Y(9);
dY(10)=4*beta*Y(9)+2*alpha*Y(11)+fco_3Ca_o*Y(7)-(3*beta+alpha+foc_3Ca_o)*Y(10);
dY(11)=3*alpha*Y(12)+3*beta*Y(10)+fco_2Ca_o*Y(6)-(2*beta+2*alpha+foc_2Ca_o)*Y(11);
dY(12)=4*alpha*ccccy+2*beta*Y(11)+fco_Ca_o*Y(5)-(beta+3*alpha+foc_Ca_o)*Y(12);  
mBK=ccccy+cccoy+ccooy+coooy+ooooy;
end

%% ODEs system
%% currents
ica=gcal*m_inf*(V-Vca);
ikdr=gk*n*(V-Vk);
ileak=gl*(V-Vl);
isk=gsk*cinf*(V-Vk);

ibk=gbk*mBK*(V-Vk);
ik =ibk + ikdr+isk;

dv= -(ica+ik+ileak)/C;
dn= (phik-n)/taun;
dc= -ff*(alfa*ica+kc*c);


if chn==1
dy=[dv dn dc dY(4) dY(5) dY(6) ]';
elseif chn==2
dy=[dv dn dc dY(4) dY(5) dY(6) dY(7) dY(8) ]';
elseif chn==3
dy=[dv dn dc dY(4) dY(5) dY(6) dY(7) dY(8) dY(9) dY(10) ]';
elseif chn==4
dy=[dv dn dc dY(4) dY(5) dY(6) dY(7) dY(8) dY(9) dY(10) dY(11) dY(12) ]';
end















