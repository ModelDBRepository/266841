function dy=model_pituitary_singleODE(t,Y,flag,struct,tag_minf)
%% A pituitary lactotroph model with BK coupled with n CaVs. 
% The I_BK current is modeled by equation (28) (BK coupled with non-inactivating  CaVs), 
% where the BK activation is modeled by equation (26)  (assuming CaV
% kinetics) or by equation (29) (assuming instantaneous activation of CaVs, tag_minf=1 in the code).
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
mBK=Y(4);
mCa=Y(5);



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
m_infV=1/(1+exp((Vm-V)/sm));
cinf=c^2/(c^2+ks^2);
alpha=m_infV/taum;
beta=(1-m_infV)/taum;
%% currents
ica=gcal*m_infV*(V-Vca);

ikdr=gk*n*(V-Vk);
ileak=gl*(V-Vl);
isk=gsk*cinf*(V-Vk);

%% BK fluxes
ca_o=struct.nca*i_ca_single/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;
ca_c=0.2;
fco_Ca_o=k1*(ca_o^n1)/(ca_o^n1+K1^n1);
foc_Ca_o=k2*(K2^n2)/(K2^n2+ca_o^n2);
foc_Ca_c=k2*(K2^n2)/(K2^n2+ca_c^n2);
k1_ca_o=fco_Ca_o;
k2_ca_o=foc_Ca_o;
k2_ca_c=foc_Ca_c;%k-C e

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

y=mBK;
%% BK 1 Ca %%%%%%%%%
if chn==1

    if tag_minf % if equal to 1 then assume that oy=my 
        m_inf_v=m_infV;
        cy=(1-m_inf_v)*y;
        oy=m_inf_v*y;
        dmBK=-k2_ca_c*cy+k1_ca_o*(m_inf_v-oy)-k2_ca_o*oy; % ode for BK activation
    else
        m_inf=mCa;
        tau_BK=(alpha+beta+k2_ca_c)/((alpha+k2_ca_c)*(k1_ca_o+k2_ca_o)+beta*k2_ca_c);
        mBK_inf=k1_ca_o*m_inf*tau_BK;
        
        dmBK=(mBK_inf-y)/tau_BK; % ode for BK activation
    end

end
%% BK 2 Ca %%chn%%%%%%%
if chn==2
    
   if tag_minf
        n_ca=2;
        m_inf_v=m_infV;

        ca2_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
        ca2_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
        ca2_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);
     
        ccy=ca2_o0_v*y;
        coy=ca2_o1_v*y;
        ooy=ca2_o2_v*y;
        dmBK=-k2_ca_c*ccy+k1_ca_o*(2*m_inf_v*(1-m_inf_v)-coy)-k2_ca_o*coy-k2_2ca_o*ooy+k1_2ca_o*(m_inf_v^2-ooy);
   else
        m_inf=mCa;
        A=beta/(2*alpha+k2_ca_c+beta);  % A_1 in the paper
        B_0=((alpha+k1_ca_o+k2_ca_o)*(1-A)+2*beta+k2_ca_c*A);   % B_2 in the paper
        B=2*beta/B_0;  % A_2 in the paper
        D=2*(1-m_inf)*m_inf*k1_ca_o/B_0;  % D_2 in the paper


        ooy=(1-B)*y-D;
        coy=(1-A)*(B*y+D);
        ccy=A*(B*y+D);
        dmBK=-k2_ca_c*ccy+k1_ca_o*(2*m_inf*(1-m_inf)-coy)-k2_ca_o*coy-k2_2ca_o*ooy+k1_2ca_o*(m_inf^2-ooy);
    end
    
    
    
end
%% BK 3 Ca %%%%%%%%%
if chn==3
   
    if tag_minf
        
        n_ca=3;
        m_inf_v=m_infV;
        ca3_o0_v=nchoosek(n_ca,0).*(1-m_inf_v).^(n_ca);
        ca3_o1_v=nchoosek(n_ca,1).*m_inf_v.*(1-m_inf_v).^(n_ca-1);
        ca3_o2_v=nchoosek(n_ca,2).*m_inf_v.^2.*(1-m_inf_v).^(n_ca-2);
        ca3_o3_v=nchoosek(n_ca,3).*m_inf_v.^3.*(1-m_inf_v).^(n_ca-3);
     
        cccy=ca3_o0_v*y;
        ccoy=ca3_o1_v*y;
        cooy=ca3_o2_v*y;
        oooy=ca3_o3_v*y;
         dmBK=-k2_ca_c*cccy+k1_ca_o*(3*m_inf_v*(1-m_inf_v)^2-ccoy)-k2_ca_o*ccoy-k2_2ca_o*cooy+k1_2ca_o*(3*m_inf_v^2*(1-m_inf_v)-cooy)-k2_3ca_o*oooy+k1_3ca_o*(m_inf_v^3-oooy);





        
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
        dmBK=-k2_ca_c*cccy+k1_ca_o*(3*m_inf*(1-m_inf)^2-ccoy)-k2_ca_o*ccoy-k2_2ca_o*cooy+k1_2ca_o*(3*m_inf^2*(1-m_inf)-cooy)-k2_3ca_o*oooy+k1_3ca_o*(m_inf^3-oooy);

    end
end


%% BK 4 Ca %%%%%%%%%
if chn==4
    
    if tag_minf
        n_ca=4;
        m_inf_v=m_infV;
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
        dmBK=-k2_ca_c*ccccy+k1_ca_o*(4* m_inf_v*(1- m_inf_v)^3-cccoy)-k2_ca_o*cccoy-k2_2ca_o*ccooy+k1_2ca_o*(6* m_inf_v^2*(1- m_inf_v)^2-ccooy)-k2_3ca_o*coooy+k1_3ca_o*(4* m_inf_v^3*(1- m_inf_v)-coooy)-k2_4ca_o*ooooy+k1_4ca_o*( m_inf_v^4-ooooy);
        
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
        dmBK=-k2_ca_c*ccccy+k1_ca_o*(4*m_inf*(1-m_inf)^3-cccoy)-k2_ca_o*cccoy-k2_2ca_o*ccooy+k1_2ca_o*(6*m_inf^2*(1-m_inf)^2-ccooy)-k2_3ca_o*coooy+k1_3ca_o*(4*m_inf^3*(1-m_inf)-coooy)-k2_4ca_o*ooooy+k1_4ca_o*(m_inf^4-ooooy);
    end
    
    
end
    


%% ODEs system
ibk=gbk*mBK*(V-Vk);
ik =ibk + ikdr+isk;

dv= -(ica+ik+ileak)/C;
dn= (phik-n)/taun;
dc= -ff*(alfa*ica+kc*c);
dmCa=(m_infV-mCa)/taum;

dy=[dv dn dc dmBK dmCa]';


