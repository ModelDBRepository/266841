% Run this script for simulating AP in a model of lactotrophs with BK coupled
% with n CaVs (here, n=nch=1) (Figure 3E). The I_BK current is
% described by equation (28) (BK coupled with non-inactivating  CaVs), 
% where the BK activation is modeled by the complete
% BK model with 2*(n+1) states (see equations (S19)-(S24)), by equation (26), 
% or by equation (29).


clc
clear all
close all



%% Pituitary model

% parameters LH

struct.C=10;  %pF
struct.gca=2;  %nS
struct.Vca=60; %mV
struct.Vm=-20; %mV
struct.sm=12;  %mV
struct.gk=3;   %nS
struct.Vk=-75; %mV
struct.vn=-5;  %mV
struct.sn=10;  %mV
struct.taun=30;%ms
struct.Vf=-20; %mV
struct.Sf=2;   %mV
struct.gl=0.2;  %nS
struct.Vl=-50; %mV
struct.fc=0.01;
struct.alfa=0.0015; %microM*fC-1
struct.Kc=0.12; %ms
struct.gsk=1.2;
struct.ks=0.4;
struct.gbk=1;

%% nCa channels and dip or ind
type='ind';
nch=1;

if nch==1
   type_color='c';
   struct.nch=1;
   struct.nca=1;
   Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1];
elseif nch==2
        type_color='g';
        if strcmp(type,'ind')
            struct.nch=2;
            struct.nca=1;
            Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1,0.1,0.1];
        else 
            struct.nch=1; 
            struct.nca=2;
            Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1];
        end
        
    elseif nch==3
        type_color='m';
        if strcmp(type,'ind')
            struct.nch=3;
            struct.nca=1;
            Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1,0.1,0.1,0.1,0.1];
            
        else
            
            struct.nch=1; 
            struct.nca=3;
            Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1];
        end
    
elseif nch==4
        type_color='r';
        if strcmp(type,'ind')
            struct.nch=4;
            struct.nca=1;
            Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1,0.1,0.1,0.1,0.1,0.1,0.1];
        else
      
            struct.nch=1; 
            struct.nca=4;
            Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1];
        end
end

    


%% ODE  system solution (equations (S19)-(S25))


tspan=0:0.05:3000;
[T, Y]=ode15s('model_pituitary_ODE',tspan, Y0,[],struct);
Y0=Y(end,:);
[T, Y]=ode15s('model_pituitary_ODE',tspan, Y0,[],struct);
color='b';

V=Y(:,1);
mBK=Y(:,2);
n=Y(:,3);
c=Y(:,4);

V_ODE=Y(:,1);
mBK_ODE=Y(:,2);
n_ODE=Y(:,3);
c_ODE=Y(:,4);

%% 1ODE solution eq (26)

Y0=[-60, 0.01, 0.1, 0.1,0.1];
if nch==1
   struct.nch=1;
   struct.nca=1;
   
elseif nch==2
        if strcmp(type,'ind')
            struct.nch=2;
            struct.nca=1;
            %Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1,0.1,0.1];
        else 
            struct.nch=1; 
            struct.nca=2;
        end
        
    elseif nch==3
        if strcmp(type,'ind')
            struct.nch=3;
            struct.nca=1;
            %Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1,0.1,0.1,0.1,0.1];
            
        else
            
            struct.nch=1; 
            struct.nca=3;
        end
    
elseif nch==4
        if strcmp(type,'ind')
            struct.nch=4;
            struct.nca=1;
            %Y0=[-60, 0.01, 0.1, 0.1, 0.1, 0.1,0.1,0.1,0.1,0.1,0.1,0.1];
        else
      
            struct.nch=1; 
            struct.nca=4;
        end
end

tspan=0:0.05:3000;
[T, Y]=ode15s('model_pituitary_singleODE',tspan, Y0,[],struct,0);
Y0=Y(end,:);
[T, Y]=ode15s('model_pituitary_singleODE',tspan, Y0,[],struct,0);



V_1ODE=Y(:,1);
T_1ODE=T;

%% 1ODE solution eq (29)
[T, Y]=ode15s('model_pituitary_singleODE',tspan, Y0,[],struct,1);
Y0=Y(end,:);
[T, Y]=ode15s('model_pituitary_singleODE',tspan, Y0,[],struct,1);

V_1ODE_simp=Y(:,1);
T_1ODE_simp=T;



mBK=Y(:,2);
n=Y(:,3);
c=Y(:,4);
C=struct.C;  %pF
% parameters
gcal=struct.gca;  %nS
Vca=struct.Vca; %mV
Vm=struct.Vm; %mV
sm=struct.sm;  %mV
gk=struct.gk; %nS
Vk=struct.Vk; %mV
Vn=struct.vn;  %mV
sn=struct.sn;
taun=struct.taun;%ms
Vf=struct.Vf; %mV
Sf=struct.Sf;   %mV
gl=struct.gl;  %nS
Vl=struct.Vl; %mV
ff=struct.fc;
alpha=struct.alfa; %microM*fC-1
kc=struct.Kc; %ms
gbk=struct.gbk;
gsk=struct.gsk;
ks=struct.ks;

%% currents
phical=1./(1+exp((Vm-V)/sm));
cinf=c.^2./(c.^2+ks^2);
ica=gcal*phical.*(V-Vca);
ibk=gbk*mBK_ODE.*(V-Vk);
ikdr=gk*n.*(V-Vk);
ileak=gl*(V-Vl);
isk=gsk*cinf.*(V-Vk);
ik =ibk + ikdr+isk;
%% BK fluxes
D_ca=250; %microm^2 s^-1
F=9.6485*10^4; %C mol^-1
conv_F=10^(-15);  % mol a M/microm^3
conv_microM=10^6; % M to microM
k_B=500; %microM^-1 * s^-1
B_tot=30; %microM

conv_V=10^(-3); %mV to V
g_ca_PQ=2.7; %pS 
g_ca_L=2; %pS 
g_ca_T=1.7; %pS 
conv_S=10^(-12); %pS to S
r_bk=30*10^(-3);% microm

i_ca_single=abs((V-Vca))*conv_V*g_ca_L*conv_S; % C/sec
ca_o=i_ca_single/(8*pi*D_ca*F*conv_F*r_bk)*exp(-r_bk/(sqrt(D_ca/k_B/B_tot)))*conv_microM;

bk_par=[  1.1093   3.3206  2.3298    0.0223    1.6012   16.5793    0.1000    0.4614];

k1_0=bk_par(1)*100;
k2_0=bk_par(2)*100;
K1=bk_par(6);
K2=bk_par(7);
n1=bk_par(3);
n2=bk_par(8);
beta0=bk_par(4);
alpha0=-bk_par(5)*beta0;
ca_c=0.2;
k1=k1_0.*exp(-alpha0*V);
k2=k2_0.*exp(-beta0*V);
foc_Ca_o=k2.*(K2^n2)./(K2^n2+ca_o.^n2);
fco_Ca_o=k1.*(ca_o.^n1)./(ca_o.^n1+K1^n1);
foc_Ca_c=k2.*(K2.^n2)./(K2.^n2+ca_c^n2);


%% figure


offset=find(V_1ODE>=V_ODE(10),1);
figure
subplot(3,1,1)
 plot(T,V_ODE,'Color', type_color)
xlim([0 1600])
ylim([-65 5])
ylabel('V [mV]')
set(gca,'FontSize',12)
subplot(3,1,2)
plot(T_1ODE(1:end-offset+1),V_1ODE(offset:end),'Color',type_color)
xlim([0 1600])
ylim([-65 5])
ylabel('V [mV]')
set(gca,'FontSize',12)
subplot(3,1,3)
plot(T,V_1ODE_simp,type_color)

%xlabel('time[ms]','FontSize',12)
%ylabel('Voltage [mV]','FontSize',12)
%legend('ODEs system','1 ODE','1 ODE simp')
%title (['gsk=' num2str(struct.gsk) ' nS, gbk=' num2str(struct.gbk) ' nS, ' num2str(nch) ' Ca channel'],'FontSize',12)
xlim([0 1600])
ylim([-65 5])
ylabel('V [mV]')
xlabel('time [ms]')
set(gca,'FontSize',12)
set(gcf,'Units','Centimeters');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);

% 
% 







