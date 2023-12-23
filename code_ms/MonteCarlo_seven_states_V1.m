function [X,t] = MonteCarlo_seven_states(Voltage_time,Voltage,X0,t_end,struct,dt)
%% Stochastic opening of Ca channel 7 states MC
%Dimension of the state space
dimX = size(X0,1);
t=0:dt:t_end;
allocated_memory=length(t);

%initialization of the states
MC0 = zeros(1,allocated_memory);
MC1 = zeros(1,allocated_memory);
MC2 = zeros(1,allocated_memory);
MC3 = zeros(1,allocated_memory);
MC4 = zeros(1,allocated_memory);
MC5 = zeros(1,allocated_memory);
MC6 = zeros(1,allocated_memory);


MC0(1) = X0(1);

s_to_ms = 1e-3;
M_to_microM = 1e-6;
alfa_0 = 3000*s_to_ms;
beta_0 = 241*s_to_ms;

gamma = 30000*s_to_ms;
delta = 11250*s_to_ms;
epsilon = 2.5e6*s_to_ms*M_to_microM;
csi = 2*s_to_ms;

struct.rate(1) = alfa_0;
struct.rate(2) = beta_0;

struct.rate(3) = gamma;
struct.rate(4) = delta;

struct.rate(5) = epsilon;
struct.rate(6) = csi;


for ind=1:allocated_memory
   
    alfa_0 = struct.rate(1);
    beta_0 = struct.rate(2);
    
    F = 96485.34;
    R = 8.314;
    T = 310; 
    if length(Voltage)~=1
    
    V = interp1(Voltage_time,Voltage,t(ind));
    else
    V=Voltage;
    end
        
    esp = (F*V*1e-3)/(R*T);
    qf = 1.15;%C
    qb = 1.94;%C
    
    % rate constants
    alfa = alfa_0*exp(qf*esp);
    beta = beta_0*exp(-qb*esp);
    

    gamma = struct.rate(3);
    delta = struct.rate(4);

   
    csi = struct.rate(6);
    
    %%%% current and calcium local %%%%%
    s_to_ms = 1e-3;
    M_to_microM = 1e-6;
   
    F=96485.34;
    R=8.314;
    T=310;
    
    
    epsilon = 2.5e6*s_to_ms*M_to_microM;

    Dca = 250; % microm*s-1;
    g_ca=2.8;
    Vca = 60;  %mV
    ica=(g_ca*(V-Vca))*1e-3; %pA
    kpiuB = 500; %microM-1*s-1
    Ca_bg = 0.1; % microM
    Btotal = 30; %microM
    r = 0.007; % microM;
    Caf=abs(ica)*1e9/(8*pi*Dca*F*r)*exp(-r/sqrt(Dca/kpiuB*Btotal))+Ca_bg;
    
    % random variables
    
    Y = rand(1);
    nReactions=12;
    Q =  zeros(1,nReactions);
    P0=MC0(ind);
    P1=MC1(ind);
    P2=MC2(ind);
    P3=MC3(ind);
    P4=MC4(ind);
    P5=MC5(ind);
    P6=MC6(ind);
    
    %% Propensity of the controller----------------------------------------
    Q(1,1) = 4*alfa*dt*P0;
    Q(1,2) = beta*dt*P1; 
    Q(1,3) = 3*alfa*dt*P1;
    Q(1,4) = 2*beta*dt*P2; 
    Q(1,5) = 2*alfa*dt*P2;
    Q(1,6) = 3*beta*dt*P3; 
    Q(1,7) = alfa*dt*P3;
    Q(1,8) = 4*beta*dt*P4;
    Q(1,9) = gamma*dt*P4;
    Q(1,10) = delta*dt*P5;
    Q(1,11) = epsilon*dt*Caf*P5;
    Q(1,12) = csi*dt*P6;
    P=[P0;P1;P2;P3;P4;P5;P6];
   
    current_state=find(P==1);
    
    if current_state==1
    
        if Y<Q(1,1) % from state 1 to state 2 MC0
            
            MC0(ind+1)=MC0(ind)-1;
            MC1(ind+1)=MC1(ind)+1;
            
        else % remain in state 1
            MC0(ind+1)=MC0(ind);
           
        end
        
    elseif current_state==2 %MC1
        
        if Y<Q(1,2)
            
            MC0(ind+1)=MC0(ind)+1;
            MC1(ind+1)=MC1(ind)-1;
            
        elseif Y<Q(1,2)+Q(1,3)
            
            MC1(ind+1)=MC1(ind)-1;
            MC2(ind+1)=MC2(ind)+1;
        else
            MC1(ind+1)=MC1(ind);
            
        end
   
     elseif current_state==3 %MC2
         if Y<Q(1,4)
            
            MC1(ind+1)=MC1(ind)+1;
            MC2(ind+1)=MC2(ind)-1;
            
        elseif Y<Q(1,4)+Q(1,5)
            
            MC2(ind+1)=MC2(ind)-1;
            MC3(ind+1)=MC3(ind)+1;
        else
            MC2(ind+1)=MC2(ind);
            
        end
         
     elseif current_state==4 %MC3
         if Y<Q(1,6)
            
            MC2(ind+1)=MC2(ind)+1;
            MC3(ind+1)=MC3(ind)-1;
            
        elseif Y<Q(1,6)+Q(1,7)
            
            MC3(ind+1)=MC3(ind)-1;
            MC4(ind+1)=MC4(ind)+1;
        else
            MC3(ind+1)=MC3(ind);
            
        end
     elseif current_state==5 %MC4
         if Y<Q(1,8)
            
            MC3(ind+1)=MC3(ind)+1;
            MC4(ind+1)=MC4(ind)-1;
            
        elseif Y<Q(1,8)+Q(1,9)
            
            MC4(ind+1)=MC4(ind)-1;
            MC5(ind+1)=MC5(ind)+1;
        else
            MC4(ind+1)=MC4(ind);
            
        end
          elseif current_state==6 %MC5
              
             if Y<Q(1,10)
            
            MC4(ind+1)=MC4(ind)+1;
            MC5(ind+1)=MC5(ind)-1;
            
        elseif Y<Q(1,10)+Q(1,11)
            
            MC5(ind+1)=MC5(ind)-1;
            MC6(ind+1)=MC6(ind)+1;
        else
            MC5(ind+1)=MC5(ind);
            
             end 
              
        elseif current_state==7
           if Y<Q(1,12)
            
            MC5(ind+1)=MC5(ind)+1;
            MC6(ind+1)=MC6(ind)-1;
            
           else
            
            MC6(ind+1)=MC6(ind);
              
           end
    end
    
    
     
end
X(1,:) = MC0(1:allocated_memory);
X(2,:) = MC1(1:allocated_memory);
X(3,:) = MC2(1:allocated_memory);
X(4,:) = MC3(1:allocated_memory);
X(5,:) = MC4(1:allocated_memory);
X(6,:) = MC5(1:allocated_memory);
X(7,:) = MC6(1:allocated_memory);
