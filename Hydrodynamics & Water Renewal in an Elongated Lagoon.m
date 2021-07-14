%% MATHEMATICAL MODELING OF HYDRODYNAMICS AND WATER QUALITY IN COASTAL REGIONS
%  Supervisor: Prof. Anastasios Stamou
%  Technical University of Munich
%  
%  Hydrodynamics and Water Renewal in an Elongated Lagoon
%  Author: Mehmet SEBER
%          github.com/mehmetseber
%          linkedin.com/in/mehmetseber/

%% Calculations
clc
clear

%  Givens for the problem
g = 9.81;                   % [m/s^2]
h = 10;                     % [m]
L = 10000;                  % [m]
n = 20;                     % [#]
B = 200;                    % [m]
tmin = 0;                   % [s]
tmax = 60000;               % [s]
T = 14400;                  % [s]
C = 50;                     % [m^0.5/s]
Zita_node = 1;              % [m]

for i = 1:n+1               % Define over the whole lagoon length
    h(i) = 10;
    B(i) = 200;
end

deltax = L/n;               % [m] distance between nodes 
c = sqrt(g*h(1));           % [m/s] celerity 
maxdeltat = deltax/c;       % [s] max. dt that can be selected 
deltat = 10;                % [s] selected deltat for smoother plots
nt = tmax/deltat;           % [#] timesteps

% Initial Conditions        % For an arbitrary matrix A ( t , i );
Q = zeros(nt+1,n+1);        %   - "t" is the # of the current row
H = zeros(nt+1,n+1);        %   - "i" is the # of the current column
A = zeros(nt+1,n+1);        % so to go rightwards in the lagoon, we
R = zeros(nt+1,n+1);        % increase i (column #) and to go forward in 
Zita = zeros(nt+1,n);       % time, we increase t (row #).

% Initial area "A", height "H" and hydraulic radius "R"
for i=1:n+1                 % Define over the whole lagoon length
    H(1,i) = h(1); 
    A(1,i) = B(1)*H(1,i);
    R(1,i) = A(1,i)/(B(1)+2*H(1,i));
end

for t=2:nt                  % Main loop for timestep
    
    time = (t-2)*deltat;    % For the first calculation, time = 0 s.
    
    % Tidal Boundary Condition Calculations
    Z_i2 = Zita_node*sin((2*pi*(time-deltat)/T)-(deltax)/(sqrt(g*H(t-1,1))*T));
    Z_i1 = Zita_node*sin((2*pi*(time-deltat))/T);
    Z_r2 = Zita(t-1,2)-Z_i2; 
    Z_r1 = Zita(t-1,1)-Z_i1;
    Z_a1 = Z_r1+(deltat/deltax)*sqrt(g*H(t-1,1))*(Z_r2-Z_r1);
    
    % Water Elevation at upstream boundary
    Zita(t,1)=Zita_node*sin((2*pi*time)/T)+Z_a1;
    
    % H, A and R calculations for the flow rate
    H(t,1) = h(1) + Zita(t,1);
    A(t,1) = B(1)*H(t,1);
    R(t,1) = A(t,1)/(B(1)+2*H(t,1));
    
    for i=2:n               % Calculations over the whole lagoon length
        
        % Water Elevations for every inside node
        Zita(t,i) = Zita(t-1,i)-(2*(deltat))/((B(i)+B(i+1))*deltax)*(Q(t-1,i+1)-Q(t-1,i));
    
        if Q(t-1,i)<0         % Negative flow case
        
            Q1 = 1/(2*deltax)*((Q(t-1,i+1))^2/A(t-1,i+1)-(Q(t-1,i-1))^2/A(t-1,i-1));
  
        else
            
            Q1 = -1/(2*deltax)*((Q(t-1,i+1))^2/A(t-1,i+1)-(Q(t-1,i-1))^2/A(t-1,i-1));
            
        end
    
        Q2 = -g*A(t-1,i)*((Zita(t,i)-Zita(t,i-1))/deltax);
        Q3 = -((g*A(t-1,i))/(C^2*R(t-1,i)))* (Q(t-1,i)/A(t-1,i))*abs((Q(t-1,i)/A(t-1,i)));
    
        % Flow Rates for every inside node
        Q(t,i) = Q(t-1,i)+Q1*deltat+Q2*deltat+Q3*deltat;
        
        % Zero Gradient Condition for the upstream boundary
        Q(t,1) = Q(t,2);
        
        % H, A and R calculations for every inside node
        H(t,i) = h(i)+(Zita(t,i)+Zita(t,i-1))/2;
        A(t,i) = B(i)*H(t,i);
        R(t,i) = A(t,i)/(B(i)+2*H(t,i));
          
    end
    
        % Downstream boundary conditions
        H(t,n+1) = h(n+1)+(Zita(t,n)+Zita(t,n))/2;
        A(t,n+1) = B(n+1)*H(t,n+1);
        R(t,n+1) = A(t,n+1)/(B(n+1)+2*H(t,n+1));
        
        % Downstream Flow Rate boundary condition
        Q(t,n+1) = Zita(t,n)*sqrt(g*B(n+1)*A(t-1,n+1));
        
end

fprintf('Calculations are done with no errors.\n'); 
%% Plots

close all
time_plot_x = [0:deltat:tmax];
 
% Water Elevation at Upstream
figure
plot(time_plot_x,Zita(:,1),'b')
title('Zita at Upstream')
xlabel('Time [s]')
ylabel('Elevation - Zita [m]')
xlim([0 tmax+deltat])
ylim([-1.1 1.1])
grid('on')
 
% Flow Rate at Upstream 
figure
plot(time_plot_x,Q(:,1))
title('Q at Upstream')
xlabel('Time [s]')
ylabel('Q [m^3/s]')
xlim([0 tmax+deltat])
grid('on')
 
% Definition for the wanted data intervals
start_1 = 43200;
End_1 = 57600;
t_interval= start_1:deltat:End_1;      
start = start_1/deltat;                   
End = End_1/deltat;                     
 
% Plot temporal variation of water elevation Zita at i = 2, 10, 15, 20
figure
hold on
plot(t_interval,(Zita(start:End,2)),'DisplayName','i=2')
plot(t_interval,(Zita(start:End,10)),'DisplayName','i=10')
plot(t_interval,(Zita(start:End,15)),'DisplayName','i=15')
plot(t_interval,(Zita(start:End,20)),'DisplayName','i=20')
legend('show')
xlim([start_1 End_1])
grid('on')
title('Zita at Locations i')
xlabel('Time [s]')
ylabel('Zita [m]')
 
 
% Plot temporal variation of flow rate Q at i = 2, 10, 15, 20
figure
hold on
plot(t_interval,(Q(start:End,2)),'DisplayName','i=2')
plot(t_interval,(Q(start:End,10)),'DisplayName','i=10')
plot(t_interval,(Q(start:End,15)),'DisplayName','i=15')
plot(t_interval,(Q(start:End,20)),'DisplayName','i=20')
legend('show')
xlim([start_1 End_1])
grid('on')
title('Q at Locations i')
xlabel('Time [s]')
ylabel('Discharge Q [m^3/s]')
 
 
% Maximum Q and Zita values
Q_max = max(Q);        
Zita_max = max(Zita);         
 
 
% Plot spatial distribution of maximum elevation 
figure
plot((0.5*deltax):deltax:L,Zita_max,'b')
xlim([0 L])
xlabel('Lenght of Lagoon [m]')
 
ylim([0 1.1])
legend('Elevation-Zita')
title('Spatial Distribution of Maximum Elevation - Zita')
ylabel('Zita [m]')
grid('on')
 
% Plot spatial distribution of maximum flow rate
figure
plot(0:deltax:L,Q_max)
xlim([0 L])
xlabel('Lenght of Lagoon [m]')
legend('Flow Rate Q')
title('Spatial Distribution of Maximum Flow Rate - Q')
ylabel('Flow Rate Q [m^3/s]')
grid('on')
