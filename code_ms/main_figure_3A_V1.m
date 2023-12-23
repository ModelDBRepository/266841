% Run this scrip for simulating AP in a neuronal model with 1:n BK-CaV
% complexes (n=1, 2 or 4) (Figure 3A). The I_BK current is modeled by
% equation (28) (BK coupled with non-inactivating  CaVs), 
% where the BK activation is modeled by equation (26).

close all
clear all
clc

% define global parameters
global gBK no_leak_na_k  gKA gKAHP gKSOR n_ca_o global_cBK_L r_bk  VmCa n_ind tau_inf tag_minf kmCa

no_leak_na_k=1; % constant leak
tag_minf=0; % no instantaneous activation of CaVs

n_ca_o=1;   % number of  synchronized CaV channels
tau_inf=1.25; % time constant for CaV

r_bk=13*10^(-3); % distance between CaV and BK channels

% current parameters

VmCa=15;
kmCa=7;
%kmCa=15;

global_cBK_L=1;
% gKAHP=0.18;
% gKA=14;

gKSOR=0;
%gKSOR=0.06;
gKAHP=0.18;
gKA=14;


% time
tspan=-30:1e-2:100;

% Initial conditions 
v_0=-67;
hNa_0=1;
mKDR_0=0.1;
mKA_0=0.1;
hKA_0=1;
mKSOR_0=0.1;
p_0=0.1;
Ca_BK_0=113;
Ca_SK_0=113;
Ca_i_0=113;
ta=3;
t0=1;
tu=[t0 t0+ta];
u_in=[-10 0];
x_0=[v_0, hNa_0, mKDR_0, mKA_0, hKA_0, mKSOR_0, p_0 , Ca_BK_0,Ca_SK_0,Ca_i_0,0,0.1];

n_sim=5; % number of tests 

for idx_sim=1:n_sim
    
    if idx_sim==1
        
        n_ind=1;    % number of CaVs
        type_line_col='-b';
        gBK=0;
    elseif idx_sim==2
        n_ind=1;    % number of CaVs
        type_line_col='-c';
        gBK=1;
    elseif idx_sim==3
        n_ind=1;    % number of CaVs
        type_line_col='-k';
        gBK=4;
    elseif idx_sim==4
        n_ind=2;    % number of CaVs
        type_line_col='-g';
        gBK=1;
    elseif idx_sim==5
        n_ind=4;    % number of CaVs
        type_line_col='-r';
        gBK=1;
    end
    %neuronal ODE model
    [T, Y]=ode15s('neuro_model_bk_cav_1_n',tspan, x_0,[],tu,u_in);
    
    %state variable v (voltage)
    v=Y(:,1);
    v_mtx(idx_sim,:)=v;
    
    f=figure (1);
    hold on
    grid on
    plot(T,v,type_line_col)
    xlim([-2 20])
end


title('Figure 3A')
%legend('iber. (g_{BK}=0)','n_{CaV}=1 g_{BK}=1','n_{CaV}=1 g_{BK}=4','n_{CaV}=2  g_{BK}=1','n_{CaV}=4 g_{BK}=1')
xlabel('time [ms]');    
ylabel('V [mV]')
    
%set(f,'Position',[10 10 100 100])
set(findall(f,'type','text'),'fontSize',8)
ylim([-80 40])
axes('Position',[.55 .55 .3 .3])
box on
hold on
grid on
idx_t1= find(T==3.75);
idx_t2=find(T==16);
plot(T(idx_t1:idx_t2),v(idx_t1:idx_t2),type_line_col);
plot(T(idx_t1:idx_t2),v_mtx(1,idx_t1:idx_t2),'-b')
plot(T(idx_t1:idx_t2),v_mtx(2,idx_t1:idx_t2),'-c')
plot(T(idx_t1:idx_t2),v_mtx(3,idx_t1:idx_t2),'-k')
plot(T(idx_t1:idx_t2),v_mtx(4,idx_t1:idx_t2),'-g')
plot(T(idx_t1:idx_t2),v_mtx(5,idx_t1:idx_t2),'-r')
ylim([-79.5 -60])
xlim([4.5 16])