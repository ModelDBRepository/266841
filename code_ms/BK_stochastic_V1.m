function [X,t] = BK_stochastic(X0,t_end,dt,struct)
%% Stochastic opening of BK channel (10 states)
%Dimension of the state space
dimX = size(X0,1);
t = 0:dt:t_end;% ms
allocated_memory = length(t);

%initialization of the states
MC0 = zeros(1,allocated_memory);
MC1 = zeros(1,allocated_memory);
MC2 = zeros(1,allocated_memory);
MC3 = zeros(1,allocated_memory);
MC4 = zeros(1,allocated_memory);
MC5 = zeros(1,allocated_memory);
MC6 = zeros(1,allocated_memory);
MC7 = zeros(1,allocated_memory);
MC8 = zeros(1,allocated_memory);
MC9 = zeros(1,allocated_memory);
MC0(1) = X0(1);
MC1(1) = X0(2);
MC2(1) = X0(3);
MC3(1) = X0(4);
MC4(1) = X0(5);
MC5(1) = X0(6);
MC6(1) = X0(7);
MC7(1) = X0(8);
MC8(1) = X0(9);
MC9(1) = X0(10);
%% reactions
s_to_ms = 1e-3;
M_to_microM = 1e-6;


KonO = 1e9*s_to_ms*M_to_microM;
KonC = 1e9*s_to_ms*M_to_microM;
KoffO = 1065*s_to_ms;
KoffC = 11917*s_to_ms;


nReactions=26;
Caloc = struct.Ca;
%rng(seed_number)
for ind=1:allocated_memory
  % random variables
    V = interp1(struct.Voltage_time,struct.Voltage,t(ind));
    Y = rand(1);
    Q = zeros(1,nReactions);
    P0 = MC0(ind);
    P1 = MC1(ind);
    P2 = MC2(ind);
    P3 = MC3(ind);
    P4 = MC4(ind);
    P5 = MC5(ind);
    P6 = MC6(ind);
    P7 = MC7(ind);
    P8 = MC8(ind);
    P9 = MC9(ind);
    
    %% time dependent variables
    
F = 96485.34;
R = 8.314;
T = 310; 
qf = 0.73;
qb = 0.58;
esp=(F*V*1e-3)/(R*T);

alfa1 = 5.5*s_to_ms*exp(qf*esp);
beta1 = 8669*s_to_ms*exp(-qb*esp);
alfa2 = 8*s_to_ms*exp(qf*esp);
beta2 = 1127*s_to_ms*exp(-qb*esp);
alfa3 = 2*s_to_ms*exp(qf*esp);
beta3 = 25.2*s_to_ms*exp(-qb*esp);
alfa4 = 884*s_to_ms*exp(qf*esp);
beta4 = 1013*s_to_ms*exp(-qb*esp);
alfa5 = 900*s_to_ms*exp(qf*esp);
beta5 = 125.7*s_to_ms*exp(-qb*esp);

    
    
    KO=1.06;
    KC=11.92;
    
    %% Propensity of the controller----------------------------------------
    Q(1,1) = 4*KonO*Caloc(ind)*dt*P0;
    Q(1,2) = beta1*dt*P0; 
    Q(1,3) = 3*KonO*Caloc(ind)*dt*P1;
    Q(1,4) = beta2*dt*P1; 
    Q(1,5) = KoffO*dt*P1;
    Q(1,6) = 2*KonO*Caloc(ind)*dt*P2;
    Q(1,7) = beta3*dt*P2;
    Q(1,8) = 2*KoffO*dt*P2;
    Q(1,9) = KonO*Caloc(ind)*dt*P3;
    Q(1,10) = beta4*dt*P3;
    Q(1,11) = 3*KoffO*dt*P3;
    Q(1,12) = beta5*dt*P4;
    Q(1,13) = 4*KoffO*dt*P4;
    Q(1,14) = 4*KonC*Caloc(ind)*dt*P5;
    Q(1,15) = alfa1*dt*P5;
    Q(1,16) = 3*KonC*Caloc(ind)*dt*P6;
    Q(1,17) = alfa2*dt*P6;
    Q(1,18) = KoffC*dt*P6;
    Q(1,19) = 2*KonC*Caloc(ind)*dt*P7;
    Q(1,20) = alfa3*dt*P7;
    Q(1,21) = 2*KoffC*dt*P7;
    Q(1,22) = KonC*Caloc(ind)*dt*P8;
    Q(1,23) = alfa4*dt*P8;
    Q(1,24) = 3*KoffC*dt*P8;
    Q(1,25) = alfa5*dt*P9;
    Q(1,26) = 4*KoffC*dt*P9;
    P=[P0;P1;P2;P3;P4;P5;P6;P7;P8;P9];
   
    current_state = find(P==1);
    
    if current_state == 1 %MC0
    
        if Y < Q(1,1) %4*KonO*KO*P0
            MC0(ind+1) = MC0(ind)-1;
            MC1(ind+1) = MC1(ind)+1;
       elseif Y < Q(1,1)+Q(1,2)% alfa1*P0
            MC0(ind+1) = MC0(ind)-1;
            MC5(ind+1) = MC5(ind)+1;
        else 
            MC0(ind+1) = MC0(ind);   
        end
        
    elseif current_state == 2 %MC1
        
        if Y < Q(1,3)%3*KonO*KO*P1
            MC1(ind+1) = MC1(ind)-1;
            MC2(ind+1) = MC2(ind)+1;
       elseif Y < Q(1,3)+Q(1,4) %beta1*P1
            MC1(ind+1) = MC1(ind)-1;
            MC6(ind+1) = MC6(ind)+1;
       elseif Y < Q(1,3)+Q(1,4)+Q(1,5) %KoffO*P1
            MC0(ind+1) = MC0(ind)+1;
            MC1(ind+1) = MC1(ind)-1;
        else
            MC1(ind+1) = MC1(ind);
        end
   
     elseif current_state == 3 %MC2
         
         if Y<Q(1,6)% 2*KonO*KO*P2
            MC2(ind+1) = MC2(ind)-1;
            MC3(ind+1) = MC3(ind)+1;
        elseif Y < Q(1,6)+Q(1,7) %gamma1*P2
            MC2(ind+1) = MC2(ind)-1;
            MC7(ind+1) = MC7(ind)+1;
         elseif Y < Q(1,6)+Q(1,7)+Q(1,8) %2*KoffO*P2
            MC2(ind+1) = MC2(ind)-1;
            MC1(ind+1) = MC1(ind)+1;
         else
            MC2(ind+1) = MC2(ind);
         end
         
     elseif current_state == 4 %MC3
         
         if Y < Q(1,9) %KonO*KO*P3
            MC4(ind+1) = MC4(ind)+1;
            MC3(ind+1) = MC3(ind)-1;
        elseif Y < Q(1,9)+Q(1,10) %delta1*P3
            MC3(ind+1) = MC3(ind)-1;
            MC8(ind+1) = MC8(ind)+1;
         elseif Y < Q(1,9)+Q(1,10)+Q(1,11) % 3*KoffO*P3
            MC3(ind+1) = MC3(ind)-1;
            MC2(ind+1) = MC2(ind)+1;
         else
            MC3(ind+1) = MC3(ind);
         end
        
     elseif current_state == 5 %MC4
         
         if Y<Q(1,12) %csi1*P4
             MC9(ind+1) = MC9(ind)+1;
             MC4(ind+1) = MC4(ind)-1;
        elseif Y < Q(1,12)+Q(1,13) % 4*KoffO*P4
             MC4(ind+1) = MC4(ind)-1;
             MC3(ind+1) = MC3(ind)+1;
         else
             MC4(ind+1) = MC4(ind);
         end
        
          elseif current_state == 6 %MC5
              
        if Y < Q(1,14)  %4*KonC*KC*P5
             MC6(ind+1) = MC6(ind)+1;
             MC5(ind+1) = MC5(ind)-1;
        elseif Y < Q(1,14)+Q(1,15) %alfa2*P5
             MC5(ind+1) = MC5(ind)-1;
             MC0(ind+1) = MC0(ind)+1;
        else
             MC5(ind+1) = MC5(ind);
        end 
              
        elseif current_state == 7
           if Y < Q(1,16) %3*KonC*KC*P6
             MC7(ind+1) = MC7(ind)+1;
             MC6(ind+1) = MC6(ind)-1;
           elseif Y < Q(1,16)+Q(1,17) %beta2*P6
             MC6(ind+1) = MC6(ind)-1;
             MC1(ind+1) = MC1(ind)+1;
           elseif Y < Q(1,16)+Q(1,17)+Q(1,18) %KoffC*P6
             MC6(ind+1) = MC6(ind)-1;
             MC5(ind+1) = MC5(ind)+1;
           else
             MC6(ind+1) = MC6(ind);
           end
           
           elseif current_state == 8 %MC7
               
           if Y < Q(1,19) %2*KonC*KC*P7
              MC8(ind+1) = MC8(ind)+1;
              MC7(ind+1) = MC7(ind)-1;
           elseif Y<Q(1,19)+Q(1,20) %gamma2*P7
              MC7(ind+1) = MC7(ind)-1;
              MC2(ind+1) = MC2(ind)+1;
           elseif Y<Q(1,19)+Q(1,20)+Q(1,21) %2*KoffC*P7
              MC7(ind+1) = MC7(ind)-1;
              MC6(ind+1) = MC6(ind)+1;
           else
              MC7(ind+1) = MC7(ind);     
           end
           
           elseif current_state==9 %MC8
               
           if Y<Q(1,22) %KonC*KC*P8
              MC9(ind+1) = MC9(ind)+1;
              MC8(ind+1) = MC8(ind)-1;
           elseif Y<Q(1,22)+Q(1,23) %delta2*P8
              MC8(ind+1) = MC8(ind)-1;
              MC3(ind+1) = MC3(ind)+1;
           elseif Y<Q(1,22)+Q(1,23)+Q(1,24) %3*KoffC*P8
              MC8(ind+1) = MC8(ind)-1;
              MC7(ind+1) = MC7(ind)+1;      
           else
              MC8(ind+1) = MC8(ind);   
           end

           elseif current_state==10 %MC9
               
           if Y<Q(1,25) %csi2*P9
              MC4(ind+1)=MC4(ind)+1;
              MC9(ind+1)=MC9(ind)-1;
           elseif Y<Q(1,25)+Q(1,26)%4*KoffC*P9
              MC9(ind+1)=MC9(ind)-1;
              MC8(ind+1)=MC8(ind)+1;
           else
              MC9(ind+1)=MC9(ind);
              
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
X(8,:) = MC7(1:allocated_memory);
X(9,:) = MC8(1:allocated_memory);
X(10,:) = MC9(1:allocated_memory);


