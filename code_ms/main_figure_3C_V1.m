% Run this script for simulating AP in human beta-cells with BK coupled
% with n L-type CaVs (n=1, 2 or 4).  The I_BK current is modeled by equation (27) 
% (with inactivating L-type CaVs), where the BK activation is modeled by equation (29)
% (assuming instantaneous activation of CaVs).

clear all 
close all

global global_cBK_T global_cBK_L global_cBK_PQ 
global gBK nmCaPQ gCaPQ VmCaPQ;
global gkatp gleak gCaL gNa  r_bk


% Currents parameters

% BK coupled with L-type
global_cBK_T=0.;
global_cBK_PQ=0;

nmCaPQ=-10;
VmCaPQ=-10;
gCaPQ=0.17;

gkatp=0.015;
gleak=0.015;
gCaL=0.14;
gNa=0.4;

% distance between CaV and BK channels
r_bk=13*10^(-3);
n_ca_o=1; %n_ca_o is the number of  synchronized CaV channels (assumed n_ca_o=1)

% time
tspan=0:1e-2:600;

%initial conditions
v_0=-49;
mkv_0=0.02;
hNa_0=0.97;
hCaL_0=0.98;
hCaT_0=1;
xERG_0=0;
yERG_0=1;
mBK_L_0=0.002;
mBK_PQ_0=0.002;
mBK_T_0=0.002;

for idx_sim=1:4
    if idx_sim==1
        gBK=0; % BK block
         global_cBK_L=0;
        x_0=[v_0, mkv_0, hNa_0, hCaL_0, hCaT_0, xERG_0, yERG_0, mBK_L_0, mBK_PQ_0, mBK_T_0];
        [T, Y]=ode15s('model_human_beta_cell_bk_cav_1_1_qss_minf',tspan, x_0,[]);
        tag_line_col='-b';
        
    elseif idx_sim==2
        gBK=1;
        global_cBK_L=.7;  
        x_0=[v_0, mkv_0, hNa_0, hCaL_0, hCaT_0, xERG_0, yERG_0, mBK_L_0, mBK_PQ_0, mBK_T_0];
        % BK coupled with n=1 L-type CaV
        [T, Y]=ode15s('model_human_beta_cell_bk_cav_1_1_qss_minf',tspan, x_0,[]);
        tag_line_col='-c';
    elseif idx_sim==3
        gBK=1;
        global_cBK_L=.5; 
        x_0=[v_0, mkv_0, hNa_0, hCaL_0, hCaT_0, xERG_0, yERG_0, mBK_L_0, mBK_PQ_0, mBK_T_0,0, 0,0];
        % BK coupled with n=2 L-type CaVs
        [T, Y]=ode15s('model_human_beta_cell_bk_cav_1_2_qss_minf',tspan, x_0,[]);
        tag_line_col='-g';
    elseif idx_sim==4
        gBK=1;
        global_cBK_L=0.4; 
        x_0=[v_0, mkv_0, hNa_0, hCaL_0, hCaT_0, xERG_0, yERG_0, mBK_L_0, mBK_PQ_0, mBK_T_0,0, 0,0,0, 0,0,0, 0,0];
        % BK coupled with n=4 L-type CaVs
        [T, Y]=ode15s('model_human_beta_cell_bk_cav_1_4_qss_minf',tspan, x_0,[]);
        tag_line_col='-r';
    end


    %state variable V
    v=Y(:,1);

    f=figure  (1);
    hold on
    grid on
    plot(T,v,tag_line_col,'Linewidth',1)
    xlabel('time [ms]')
    ylabel('V [mV]')
    %set(f,'Position',[10 10 125 125])
    ylim([-75 2])
    xlim([0 550])
    title('Figure 3C')
end

