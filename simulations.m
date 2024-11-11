
%% CSTR Simulation
num_hours = 8/5;
rng(2)
 %% Define parameters
% Process parameters
tic
p.k0  = 7.2E10;             % 1/min, Arrhenius constant
p.A   = 0.1666;             % m2, Cross sectional area
u.E_R  = 8750;         % K, Activation energy divided by Universal Gas Constant
p.H_rxn = -5E4;             % J/mol, Heat of reaction
p.rho_Cp = 239;             % J/(L.K), Reactants volumetric heat capacity 
p.rho_c_Cpc = 4175;         % J/(L.K), Coolant volumetric heat capicity
p.Vc = 10;                  % L, Coolant volume
u.UAc = 5E4;                % J/(min.K), Heat transfer coefficient
u.cV_QC = 150;              % L/(min.m1/2), Coolant Valve coefficient
u.cV_Q = 150;               % L/(min.m1/2), Coolant Valve coefficient
u.T_SP = 309;               % K, Reactor temperature setpoint
p.Q_SP = 100;               % L/min, Outlet flowrate setpoint
p.h_SP = @(t) 1.6;          % m, Tank level setpoint
p.QC_SP = 15;               % L/min, Coolant flowrate setpoint
p.delta_PCU_PCD = 25;       % psi, coolant valve pressure differential
p.detla_PU_PD = 50;         % psi, outlet valve pressure differential
p.T = 1;                    % s, measurement interval
p.T_end = num_hours*3600;   % s, run time
p.CA_var = 5.56E-6;         % mol/L, noise on reactor concentration measurement
p.T_var = 4.44E-1;          % K, noise on reactor temperature measurement
p.TC_var = 4.44E-1;         % K, noise on coolant temperature measurement
p.h_var = 1.6E-3;           % m, noise on level measurement
p.Q_var = 4.44E-1;          % L/min, noise on outlet flowrate measurement
p.QC_var = 1.00E-2;         % L/min, noise on coolant flowrate measurement
p.QF_var = 4.44E-1;         % L/min, noise on inlet flowrate measurement
p.CAF_var = 5.56E-6;        % mol/L, noise on inlet concentration measurement
p.TF_var = 4.44E-1;         % K, noise on inlet temperature measurement
p.TCF_var = 4.44E-1;        % K, noise on coolant inlet temperature measurement

%% Valve characteristics
p.tau = 2;                  % Valve time constant
p.K = 1;                    % Valve gain

%% Control parameters
% Outlet valve
p.Kp_Q_valve = -0.01;          % Proportional controller gain for outlet valve
p.Ki_Q_valve = 0.1;               % Integral controller gain for outlet valve
p.Kp_Q_SP = 1;               % Proportional controller gain for flowrate setpoint
p.Ki_Q_SP = 0.1;                  % Integral controller gain for flowrate setpoint

% Coolant valve
p.Kp_QC_valve = -0.01;         % Proportional controller gain for coolant valve
p.Ki_QC_valve = 0.1;            % Integral controller gain for coolant valve
p.Kp_QC_SP = 1;              % Proportional controller gain for coolant flowrate setpoint
p.Ki_QC_SP = 0.1;               % Integral controller gain for coolant flowrate setpoint


% Valve stiction
MV_low_bound = 0;
MV_high_bound = 1;
% model parameters
S = 5; % deadband plus stick band parameter
J = 0; % slip jump parameter
currentTimeStamp = 1;
%% Define inputs and manipulated variables
% Inputs
u.QF  = @(t) 100;               % L/min, Feed flowrate
u.CAF = @(t) 1.0;               % mol/L, Reactant A feed concentration
u.QC  = @(t) 15 + 0*t;          % L/min, Coolant flowrate
u.TF  = @(t) 320;               % K, Reactor feed temperature
u.TCF = @(t) 300;               % K, Coolant feed temperature

% Measurement biases
u.QC_meas = 0;
u.T_meas = 0;

% Valve stiction
u.QC_valve_stiction = 0;

% Level control loop
u.Q_valve_PI = 0.5;             % ~, Fraction outlet valve opening
Q_valve_PI = [];
Q_valve_PI(1) = u.Q_valve_PI;   % ~, Fraction valve opening vector (keeps track of all fraction valve openings for Q_valve)
u.Q_SP = p.Q_SP;                % L/min, Outlet flowrate setpoint
Q_SP = [];
Q_SP(1) = u.Q_SP;

% Temperature control loop
u.QC_valve_PI = 0.5;            % ~, Fraction outlet valve opening
QC_valve_PI = [];               
QC_valve_PI(1) = u.QC_valve_PI; % ~, Fraction valve opening vector (keeps track of all fraction valve openings for QC_valve)
u.QC_SP = p.QC_SP;              % L/min, Outlet flowrate setpoint
QC_SP = [];
QC_SP(1) = u.QC_SP;
%% Generate operating conditions
% Modes:
% 1: Normal operation: Variable: None, Nominal value: +0
% 2: Catalyst deactivation: Variable: E_R, Nominal value: +3 
% 3: Heat exchanger fouling: Variable: UAc, Nominal value: -125
% 4: Dead coolant flow measurement: Variable: QC_meas, Nominal value: +0
% 5: Bias in temperature measurement: Variable: T_meas, Nominal value: +4
% 6: Bias in temperature measurement: Variable: T_meas, Nominal value: -4
% 7: Coolant valve stiction: Variable: QC_valve_stiction, Nominal value: +0
% 8: Coolant valve stiction: Variable: QC_valve_stiction, Nominal value: +0
% 9: Ramp change in CAF: Variable: CAF, Nominal value: +6E-4
% 10: Ramp change in CAF: Variable: CAF, Nominal value: -6E-4
% 11: Ramp change in TF: Variable: TF, Nominal value: +0.1
% 12: Ramp change in TF: Variable: TF, Nominal value: -0.1
% 13: Ramp change in TCF: Variable: TCF, Nominal value: +0.1
% 14: Ramp change in TCF: Variable: TCF, Nominal value: -0.1
% 15: Step change in PCU: Variable: cV_QC, Nominal value: +50
% 16: Step change in PCU: Variable: cV_QC, Nominal value: -50
% 17: Step change in PD: Variable: cV_Q, Nominal value: +100
% 18: Step change in PD: Variable: cV_Q, Nominal value: -100
% 19: Setpoint change in T: Variable: T_SP, Nominal value: +3
% 20: Setpoint change in T: Variable: T_SP, Nominal value: -3
% 21: Step change in QF: Variable: QF, Nominal value: +10
% 22: Step change in QF: Variable: QF, Nominal value: -10
% 23: Damped oscillations in the feed flowrate: Variable: QF, Nominal value: +10
% 24: Autoregressive disturbance in the feed flowrate: Variable: QF, Nominal value: +0
% 25: High frequency oscillations in the feed flowrate: Variable: QF, Nominal value: +10
% 26: Intermediate frequency oscillations in the feed flowrate: Variable: QF, Nominal value: +10
% 27: Intermediate frequency oscillations in the feed flowrate: Variable: QF, Nominal value: +10
% 28: High frequency oscillations in the feed flowrate: Variable: QF, Nominal value: +10
% Probabilities of faults and NOC2
% N2 = 0.001;
% a1 = 0.001;
% a2 = 0.001;
% a3 = 0.001;
% a4 = 0.001;
% a5 = 0.001;
% M2 = N2;
% M3 = a1;
% M4 = a2;
% M5 = a3;
% M6 = a4;
% M7 = a5;
% M8 = a1*a2;
% M9 = a1*a3;
% M10 = a4*a2;
% M11 = a4*a3;
% M12 = a5*a2;
% M13 = a5*a3;
% stay1 = 1;
% stay2 = 1;
% stay3 = 1;
% stay4 = 1;
% stay5 = 1;
% stay6 = 1;
% stay7 = 1;
% stay8 = 1;
% stay9 = 1;
% stay10 = 1;
% stay11 = 1;
% stay12 = 1;
% stay13 = 1;
% % Mode transition matrix
% A = [stay1 0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01    0.01   0.01   0.01;...
%      M2    stay2 M3    M4    M5    M6    M7    M8    M9    M10     M11    M12    M13;...
%      M3    M3    stay3 M3    M3    M3    M3    M4    M5    M3      M3     M3     M3;...  
%      M4    M4    M4    stay4 M4    M4    M4    M3    M4    M6      M4     M7     M4;...
%      M5    M5    M5    M5    stay5 M5    M5    M5    M3    M5      M6     M5     M7;...
%      M6    M6    M6    M6    M6    stay6 M6    M6    M6    M4      M5     M6     M6;...
%      M7    M7    M7    M7    M7    M7    stay7 M7    M7    M7      M7     M4     M5;...
%      M8    M8    M4    M3    M8    M8    M8    stay8 M8    M8      M8     M8     M8;...
%      M9    M9    M5    M9    M3    M9    M9    M9    stay9 M9      M9     M9     M9;...
%      M10   M10   M10   M6    M10   M4    M10   M10   M10   stay10  M10    M10    M10;...
%      M11   M11   M11   M11   M6    M5    M11   M11   M11   M11     stay11 M11    M11;...
%      M12   M12   M12   M7    M12   M12   M4    M12   M12   M12     M12    stay12 M12;...
%      M13   M13   M13   M13   M7    M13   M5    M13   M13   M13     M13    M13    stay13];
% stay1 = 1 - sum(A(A(:,1)~=stay1,1));
% stay2 = 1 - sum(A(A(:,2)~=stay2,2));
% stay3 = 1 - sum(A(A(:,3)~=stay3,3));
% stay4 = 1 - sum(A(A(:,4)~=stay4,4));
% stay5 = 1 - sum(A(A(:,5)~=stay5,5));
% stay6 = 1 - sum(A(A(:,6)~=stay6,6));
% stay7 = 1 - sum(A(A(:,7)~=stay7,7));
% stay8 = 1 - sum(A(A(:,8)~=stay8,8));
% stay9 = 1 - sum(A(A(:,9)~=stay9,9));
% stay10 = 1 - sum(A(A(:,10)~=stay10,10));
% stay11 = 1 - sum(A(A(:,11)~=stay11,11));
% stay12 = 1 - sum(A(A(:,12)~=stay12,12));
% stay13 = 1 - sum(A(A(:,13)~=stay13,13));
% A = [stay1 0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01    0.01   0.01   0.01;...
%      M2    stay2 M3    M4    M5    M6    M7    M8    M9    M10     M11    M12    M13;...
%      M3    M3    stay3 M3    M3    M3    M3    M4    M5    M3      M3     M3     M3;...  
%      M4    M4    M4    stay4 M4    M4    M4    M3    M4    M6      M4     M7     M4;...
%      M5    M5    M5    M5    stay5 M5    M5    M5    M3    M5      M6     M5     M7;...
%      M6    M6    M6    M6    M6    stay6 M6    M6    M6    M4      M5     M6     M6;...
%      M7    M7    M7    M7    M7    M7    stay7 M7    M7    M7      M7     M4     M5;...
%      M8    M8    M4    M3    M8    M8    M8    stay8 M8    M8      M8     M8     M8;...
%      M9    M9    M5    M9    M3    M9    M9    M9    stay9 M9      M9     M9     M9;...
%      M10   M10   M10   M6    M10   M4    M10   M10   M10   stay10  M10    M10    M10;...
%      M11   M11   M11   M11   M6    M5    M11   M11   M11   M11     stay11 M11    M11;...
%      M12   M12   M12   M7    M12   M12   M4    M12   M12   M12     M12    stay12 M12;...
%      M13   M13   M13   M13   M7    M13   M5    M13   M13   M13     M13    M13    stay13];
% C = cumsum(A,1);
% mode = 1;
% modes= [];
% r = [];
% r(1) = 1;
% modes(1) = mode;
mode_fields = ["E_R", "UAc", "QC_meas", "T_meas", "QC_valve_stiction",...
                "CAF", "TF", "TCF", "cV_QC", "cV_Q", "T_SP", "QF"];
freq = 1440; % Every operating mode is 120 minutes

nom_values = [ +0; +3; -125; +0; +4; -4; +0; +0; ...
               +6E-4; -6E-4; +0.1; -0.1; +0.1; ... 
               -0.1; +50; -50; +100; -100;...
                +3; -3; +10; -10; +10; +0; ... 
               +10; +10; +10; +10];
NOC_values = [u.E_R u.UAc u.QC_meas u.T_meas u.QC_valve_stiction...
              u.CAF(1) u.TF(1) u.TCF(1) u.cV_QC u.cV_Q u.T_SP u.QF(1)];
rng(1)
uu = zeros(p.T_end,length(NOC_values));
num_modes = 1:13;
weight_modes = ones(size(num_modes)); weight_modes(1) = 3;
num_windows = floor(p.T_end/freq);
windows = zeros(1,num_windows);
windows(1) = 1;
windows(2:end) = randsample(num_modes,num_windows-1,true,weight_modes);

windows = [1 5 5 1];
modes = windows;



for i = 1:num_windows
%     % Current mode
%     r(i) = rand;
%     mode = find(C(:,mode) >= r(i), 1);
    

    

%     if modes(i) ~= modes(i-1)
%                rand_mag = (rand + 0.25);
%                ramp = (i-1);
%     end
    mode = modes(i);
    switch mode
        case {1} % NOC1
        
        case {2} % Setpoint change, NOC2
            var = 11;
            rand_mag = (rand + 0.25);
            uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(19);
            
        case {3} % F1: Heat exchanger fouling
            var = 2;
            rand_mag = (rand + 0.25);
            ramp = 1:1024;
            uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(3).*(ramp)/60;

        case {4} % F2: Dead coolant flowrate measurement
            var = 3;
            rand_mag = (rand + 0.25);
            uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;
                
        case {5} % F3: Coolant valve stiction
            var = 5;
            rand_mag = (rand + 0.25);
            uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;

        case {6} % F4: Step change in pressure in outlet line
             var = 10;
             rand_mag = (rand + 0.25);
             uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(17);

        case {7} % F5: Auto-regressive disturbance in feed flowrate
             var = 12;
             rand_mag = (rand + 0.25);
             ramp = 1:1024;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) =  0.8*uu(i-1,var) + randn;
  
        case {8} % (F1 + F2): Heat exchanger fouling + Dead coolant flowrate measurement
             var = 2;
             rand_mag = (rand + 0.25);
             ramp = 1:1024;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(3).*(ramp)/60;
             var = 3;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;

        case {9} % (F1 + F3): Heat exchanger fouling + Coolant valve stiction
             var = 2;
             rand_mag = (rand + 0.25);
             ramp = 1:1024;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(3).*(ramp)/60;
             var = 5;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;

        case {10} % (F4 + F2):
            % Step change in pressure in outlet line + Dead coolant flowrate measurement
             var = 10;
             rand_mag = (rand + 0.25);
             uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(17);
             var = 3;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;

        case {11} % (F4 + F3): Step change in pressure in outlet line + Coolant valve stiction
             var = 10;
             rand_mag = (rand + 0.25);
             uu(((i-1)*freq + 1):((i-1)*freq + 1024), var) = rand_mag*nom_values(17);
             var = 5;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;

        case {12} % (F5 + F2): Auto-regressive disturbance in feed flowrate + Dead coolant flowrate measurement
             var = 12;
             for j = 1 : 1024
                uu((i-1)*freq + j, var) =  0.8*uu((i-1)*freq + j - 1,var) + randn;
             end
             var = 3;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;


        case {13} % (F5 + F3): Auto-regressive disturbance in feed flowrate + Coolant valve stiction
             var = 12;
             for j = 1 : 1024
                uu((i-1)*freq + j, var) =  0.8*uu((i-1)*freq + j - 1,var) + randn;
             end
             var = 5;
             uu(((i-1)*freq + 1):((i-1)*freq + 1024),var) = 1;
    end
end
uu = uu + NOC_values;
for i = 1 : length(mode_fields)
    u.(mode_fields(i)) = griddedInterpolant(linspace(1, p.T_end, length(uu)), uu(:,i));
end
%% Obtain initial steady state
% Initial process conditions
p.state_fields = {'CA', 'T', 'TC', 'h' ,'Q_valve' ,'QC_valve'};     % Field names for each state
x0.CA = 0.037;                                                      % mol/L, initial tank A concentration
x0.T = 320;                                                         % K, initial tank temperature
x0.TC = 345.44;                                                     % K, initial coolant temperature
x0.h = 1.6;                                                         % mol, initial tank liquid height
x0.Q_valve = Q_valve_PI;
x0.QC_valve = QC_valve_PI;
x0_vec = struct2vec(x0, p.state_fields);
x0 = fsolve(@(x_vec) SystemODEs(0, x_vec, u, p), x0_vec);
%% Initial conditions
x0_vec = x0;
p.h_SP = @(t) x0_vec(4) + 0*(t > 10);
% Initial control conditions
int_Err_h = 0;
int_Err_Q_SP = 0;
int_Err_T = 0;
int_Err_QC_SP = 0;

%% Simulation
sol = ode23s(@(t, x_vec) SystemODEs(t, x_vec, u, p), [0 p.T], x0_vec);
y = [];
T = p.T;
y(1,:) = Measurement(T, sol, u, p);

% Update level control loop
int_Err_h  = int_Err_h  + (y(1,4) - p.h_SP(T)) * p.T;
u.Q_SP = p.Q_SP + p.Kp_Q_SP*(y(1,4) - p.h_SP(T)) + p.Kp_Q_SP*p.Ki_Q_SP*int_Err_h;
int_Err_Q_SP = int_Err_Q_SP + (y(1,5) - u.Q_SP);
u.Q_valve_PI  = max(0, min(1, Q_valve_PI(1) + p.Kp_Q_valve*(y(1,5) - u.Q_SP) + p.Kp_Q_valve*p.Ki_Q_valve*int_Err_Q_SP));
% Update temperature control loop
int_Err_T  = int_Err_T  + (y(1,2) - u.T_SP(T)) * p.T;
u.QC_SP = p.QC_SP + p.Kp_QC_SP*(y(1,2) - u.T_SP(T)) + p.Kp_QC_SP*p.Ki_QC_SP*int_Err_T;
int_Err_QC_SP = int_Err_QC_SP + (y(1,6) - u.QC_SP);
u.QC_valve_PI  = max(0, min(1, QC_valve_PI(1) + p.Kp_QC_valve*(y(1,6) - u.QC_SP) + p.Kp_QC_valve*p.Ki_QC_valve*int_Err_QC_SP));
for i = 2:(p.T_end/p.T)
    sol = odextend(sol, @(t, x_vec) SystemODEs(t, x_vec, u, p), i*p.T);
    T = i*p.T;
    y(i,:) = Measurement(T, sol, u, p);
    if rem(i,43*freq) == 0
        save("CSTR_re_run_pre_data.mat","y","modes");
    end

    % Dead coolant flowrate measurement
    if u.QC_meas(T) ~= 0
       y(i,6) = y(i - 1,6);
    end
       

    % Anti-reset windup
    if (y(i,11) == 0 || y(i,11) == 1)
        int_Err_h = 0;
        int_Err_Q_SP = 0;
    end

    if (y(i,13) == 0 || y(i,13) == 1)
        int_Err_T = 0;
        int_Err_QC_SP =0;
    end
    % Update level control loop
    int_Err_h  = int_Err_h  + (y(i,4) - p.h_SP(T)) * p.T;
    u.Q_SP = p.Q_SP + p.Kp_Q_SP*(y(i,4) - p.h_SP(T)) + p.Kp_Q_SP*p.Ki_Q_SP*int_Err_h;
    int_Err_Q_SP = int_Err_Q_SP + (y(i,5) - u.Q_SP);
    u.Q_valve_PI  = max(0, min(1, Q_valve_PI(1) + p.Kp_Q_valve*(y(i,5) - u.Q_SP) + p.Kp_Q_valve*p.Ki_Q_valve*int_Err_Q_SP));
    Q_valve_PI(i) = u.Q_valve_PI;
    Q_SP(i) = u.Q_SP;
   
    % Update temperature control loop
    int_Err_T  = int_Err_T  + (y(i,2) - u.T_SP(T)) * p.T;
    u.QC_SP = p.QC_SP + p.Kp_QC_SP*(y(i,2) - u.T_SP(T)) + p.Kp_QC_SP*p.Ki_QC_SP*int_Err_T;
    int_Err_QC_SP = int_Err_QC_SP + (y(i,6) - u.QC_SP);
    u.QC_valve_PI  = max(0, min(1, QC_valve_PI(1) + p.Kp_QC_valve*(y(i,6) - u.QC_SP) + p.Kp_QC_valve*p.Ki_QC_valve*int_Err_QC_SP));
    QC_valve_PI(i) = u.QC_valve_PI;
    if u.QC_valve_stiction(T) ~= 0
       %% stiction calculation starts
       %% initialize actions in MV units
       if currentTimeStamp == 1
          crntAction = QC_valve_PI(i);
          nxtAction = crntAction;
       elseif currentTimeStamp > 1
              crntAction = QC_valve_PI(i - 1);
              nxtAction = QC_valve_PI(i);
       end
       %% stiction model applied to valve output
       if currentTimeStamp == 1
          xss = 0; % initialize memory variable for output signal when valve becomes stuck
          % previous output as a % of the final element range
          MV_output_previous = ( (crntAction - MV_low_bound)/ (MV_high_bound - MV_low_bound) )*100; %crntAction;
          MV_incoming_present = crntAction; % initialize present MV signal
          MV_incoming_previous = crntAction;% initialize previous MV signal
          vnew_prev = 0; % initialize previous local gradient of control signal
          I = 0; % initialize flag for valve stuck during trajectory
       else
          MV_output_previous = MV_output; % previous output as a % of the final element range
          MV_incoming_present = nxtAction;% present MV signal
          MV_incoming_previous = crntAction;% previous MV signal
          vnew_prev = vnew; % previous local gradient of control signal
        end
        %% stiction calculation, MV_output is a % of MV range
        [MV_output,I,xss,vnew,MV_Percentage_present] = stictionModel(J,S,xss,...
        MV_output_previous,I,...
        MV_incoming_present,...
        MV_incoming_previous,...
        MV_low_bound,MV_high_bound,...
        vnew_prev);
        controlLawSignal(currentTimeStamp) = QC_valve_PI(i);
        % convert MV_output at current timestamp to actual MV value
        crntAction = (MV_output/100)*(MV_high_bound - MV_low_bound) + MV_low_bound; % output after stiction model at current timestamp
        % apply stiction model output as current action
        finalElementAdjustment(currentTimeStamp) = crntAction;
        QC_valve_PI(i) = crntAction;
        u.QC_valve_PI = crntAction;
        currentTimeStamp = currentTimeStamp + 1;
    else
        currentTimeStamp = 1;
        QC_valve_PI(i) = u.QC_valve_PI;
    end
    QC_SP(i) = u.QC_SP;
end

%% Plot results
% plot(1:p.T_end,u.QF(1:p.T_end))
% title("Inlet flowrate")
% xlabel("Time [s]")
% ylabel("Flowrate [L/min]")
% ylim([0 200])
% plot(p.T:p.T:p.T_end,y(:,6),'k.')
% title("Coolant flowrate")
% xlabel("Time [s]")
% ylabel("Flowrate [L/min]")
% 
% 
% plot(sol.x,sol.y(4,:),'-',p.T:p.T:p.T_end,y(:,4),'.',1:p.T_end,p.h_SP(1:p.T_end),'k--')
% title("Level response")
% xlabel("Time [s]")
% ylabel("Level [m]")
% legend("Actual","Measurement","Setpoint")
% 
% 
% legend("Position", [0.67031,0.15346,0.22321,0.11905])
% plot(p.T:p.T:p.T_end,Q_SP,'b-',p.T:p.T:p.T_end,y(:,12),'k.')
% title("Outlet valve setpoint")
% xlabel("Time [s]")
% ylabel("Flowrate [L/min]")
% legend("Setpoint","Measurement")
% plot(sol.x,sol.y(5,:),'-',p.T:p.T:p.T_end,Q_valve_PI,'k.',p.T:p.T:p.T_end,y(:,11),'.')
% title("Outlet valve opening")
% xlabel("Time [s]")
% ylabel("Fraction valve opening [~]")
% legend("Delayed response","Control output","Measurement")
% plot(sol.x,sol.y(2,:),'-',p.T:p.T:p.T_end,y(:,2),'.',1:p.T_end,u.T_SP(1:p.T_end),'k--')
% title("Temperature response")
% xlabel("Time [s]")
% ylabel("Temperature [K]")
% legend("Actual","Measurement","Setpoint")
% 
% legend("Position", [0.67374,0.20832,0.22321,0.11905])
% plot(p.T:p.T:p.T_end,QC_SP,'b-',p.T:p.T:p.T_end,y(:,14),'k.')
% title("Coolant valve setpoint")
% xlabel("Time [s]")
% ylabel("Flowrate [L/min]")
% legend("Setpoint","Measurement")
% 
% plot(sol.x,sol.y(6,:),'-',p.T:p.T:p.T_end,QC_valve_PI,'k.',p.T:p.T:p.T_end,y(:,13),'.')
% title("Coolant valve opening")
% xlabel("Time [s]")
% ylabel("Fraction valve opening [~]")
% legend("Delayed response","Control output","Measurement")

%% Cluster modes
std_data = (y-mean(y))./std(y);
% % Perform PCA using built-in function
[coeff, score, latent, ~, explained] = pca(std_data);
C_T = modes(1:p.T:end)';
% run = {"7 days"}; 

% saving data
data.score = score;
data.explained = explained;
data.C_T = C_T;
% data.run = run;
% data.A = A;
data.time = toc;

save("CSTR_re_run.mat","data");
% K = 4;
% tol = 1e-2;
% J_best = -1000;
% C_best = [];
% pi_best = [];
% mu_best = [];
% sigma_best = [];
% A_kl_best = [];
% for PC = 1:14
%     for i = 1
%         rng(i)
%         [C_TCGMM, pi_TCGMM, mu_TCGMM, sigma_TCGMM, A_kl ] = TCGMM(score(:,1:PC), K, tol);
%         A_kl = A_kl./(sum(A_kl,1));
%           J_TCGMM = sum(diag(confusionmat(C_T, C_TCGMM)));
% %         J_TCGMM = sum_cluster_transitions(C_TCGMM);
% %         J_TCGMM = similarity(score(:,1:PC),C_TCGMM);
%         if J_TCGMM > J_best
%             J_TCGMM = J_best;
%             C_best = C_TCGMM;
%             pi_best = pi_TCGMM;
%             mu_best = mu_TCGMM;
%             sigma_best = sigma_TCGMM;
%             A_kl_best = A_kl;
%         end
%    
%     end
%     PC
% end

%%

% C_TCGMM = C_best;
% [C_TCGMM, A_kl_best] = ac_label(C_T, C_TCGMM,A_kl_best);
% subplot(1,2,1)
% gscatter(score(:,1),score(:,2),C_T,'brgc','o')
% axis equal
% title("Ground truth")
% xlabel("PC1")
% ylabel("PC2")
% subplot(1,2,2)
% gscatter(score(:,1),score(:,2),C_TCGMM,'brgc','o')
% axis equal
% title("TCGMM")
% xlabel("PC1")
% ylabel("PC2")
% figure
% t = 1:p.T_end;
% t_constant = ones(size(t))
% plot(t(C_TCGMM==1),score(C_TCGMM==1,1),'b*',t(C_TCGMM==2),score(C_TCGMM==2,1),'r*',t(C_TCGMM==3),score(C_TCGMM==3,1),'g*',t(C_TCGMM==4),score(C_TCGMM==4,1),'c*');
% plot(t(C_T==1),score(C_T==1,1)-5,'b*',t(C_T==2),score(C_T==2,1)-5,'r*',t(C_T==3),score(C_T==3,1)-5,'g*',t(C_T==4),score(C_T==4,1)-5,'c*');
% figure
% cm = confusionchart(C_T, C_TCGMM);
% cm.RowSummary = 'row-normalized';
% cm.Title = '% Accuracy';
% cm.FontSize = 18;
% perc_diff = sum(sum(abs(A_kl_best-A)*100))



%% Functions
function dxdt = SystemODEs(t, x_vec, u, p)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = vec2struct(x_vec, p.state_fields);
v = CalculateIntermediates(t, x, u, p);

% Calculate state derivatives as structure
ddt.CA =  (v.nA_in - v.nA_reac - v.nA_out)./(p.A.*x.h);
ddt.T = v.T_in + v.HRX - v.T_out + v.Heat_transfer_A;
ddt.TC = v.TC_in - v.TC_out + v.Heat_transfer_C;
ddt.h = (u.QF(t) - v.Q)./p.A;
ddt.Q_valve = ( u.Q_valve_PI - x.Q_valve) / p.tau;
ddt.QC_valve = ( u.QC_valve_PI - x.QC_valve) / p.tau;

% Map state derivative structure to vector
dxdt = struct2vec(ddt, p.state_fields);
end

function v = CalculateIntermediates(t, x, u, p)
% Calculate intermediate process variables 
% (i.e. all variables which are neither exogeneous inputs 
%       nor state variables)
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Calculate all flowrates into / out of CSTR
v.Q = x.Q_valve*u.cV_Q(t).*sqrt(x.h);                                     % L/min, flowrate out of CSTR
v.QC = x.QC_valve*u.cV_QC(t).*(2/150)*u.QC(t);
v.nA_in = u.QF(t).*u.CAF(t);                                         % mol/min, molar flowrate of A into CSTR
v.nA_out = v.Q.*x.CA;                                                % mol/min, molar flowrate of A exiting CSTR

% Calculate the amount of reactant A reacted
v.nA_reac = p.k0.*exp(-u.E_R(t)./x.T).*x.CA.*p.A.*x.h;                  % mol/min, moles of A that reacted

% Calculate energy balance-related variables
v.HRX = (v.nA_reac.*(-p.H_rxn))./(p.rho_Cp.*p.A.*x.h);               % K/min, Heat gained by reaction
v.Heat_transfer_A = (u.UAc(t).*(x.TC - x.T)) ./ (p.rho_Cp.*p.A.*x.h);   % K/min, Temperature decrease of A due to cooling of C
v.Heat_transfer_C = (u.UAc(t).*(x.T - x.TC)) ./ (p.rho_c_Cpc.*p.Vc);    % K/min, Temperature increase of C due to cooling A
v.T_in = u.QF(t).*u.TF(t)./(p.A.*x.h);                               % K/min, Heat of A entering the system
v.T_out = v.Q.*x.T./(p.A.*x.h);                                      % K/min, Heat of A exiting the system
v.TC_in = v.QC.*u.TCF(t)./p.Vc;                                   % K/min, Heat of C entering the jacket
v.TC_out = v.QC.*x.TC./p.Vc;                                      % K/min, Heat of C exiting the jacket
end

function y = Measurement(t, sol, u, p)
    x.CA = sol.y(1,end);
    x.T = sol.y(2,end);
    x.TC = sol.y(3,end);
    x.h = sol.y(4,end);
    x.Q_valve = sol.y(5,end);
    x.QC_valve = sol.y(6,end);
    v = CalculateIntermediates(t, x, u, p);
    y(1) = sol.y(1,end) + p.CA_var*randn;
    y(2) = sol.y(2,end) + p.T_var*randn + u.T_meas(t);
    y(3) = sol.y(3,end) + p.TC_var*randn;
    y(4) = sol.y(4,end) + p.h_var*randn;
    y(5) = v.Q(end) + p.Q_var*randn;
    y(6) = v.QC(end) + p.QC_var*randn;
    y(7) = u.QF(t) + p.QF_var*randn;
    y(8) = u.CAF(t) + p.CAF_var*randn;
    y(9) = u.TF(t) + p.TF_var *randn;
    y(10) = u.TCF(t) + p.TCF_var*randn;
    y(11) = x.Q_valve;
    y(12) = u.Q_SP;
    y(13) = x.QC_valve;
    y(14) = u.QC_SP;
end

function[MV_output,I,xss,vnew,MV_Percentage_present] = stictionModel(J,S,...
xss,MV_output_previous,...
I,MV_incoming_present,...
MV_incoming_previous,...
MV_low_bound,...
MV_high_bound,vnew_prev)
%% saturated condition and scaling of present MV input signal
if MV_incoming_present <= MV_low_bound
MV_Percentage_present = 0;
MV_output = 0;
elseif MV_incoming_present >= MV_high_bound
MV_Percentage_present = 100;
MV_output = 100;
elseif ( MV_incoming_present < MV_high_bound ) && ( MV_low_bound < MV_incoming_present )
MV_Percentage_present = ( (MV_incoming_present - MV_low_bound)/(MV_high_bound - MV_low_bound) )*100;
end
%% scaling of previous MV input signal
if MV_incoming_previous <= MV_low_bound
MV_Percentage_previous = 0;
elseif MV_incoming_previous >= MV_high_bound
MV_Percentage_previous = 100;
elseif ( MV_incoming_previous < MV_high_bound ) && ( MV_low_bound < MV_incoming_previous ) 
MV_Percentage_previous = ( (MV_incoming_previous - MV_low_bound)/(MV_high_bound - MV_low_bound) )*100;
end
vnew = (MV_Percentage_present - MV_Percentage_previous)/(1); % gradient of incoming control signal (consistently i.t.o. MDP time)
%%
if MV_Percentage_present > 0 && MV_Percentage_present < 100
vold = vnew_prev;
vnew = (MV_Percentage_present - MV_Percentage_previous)/(1); % gradient of incoming control signal (consistently i.t.o. MDP time)
if sign(vnew) == sign(vold)
if I ~= 1 % if not stuck during valve transient behaviour
% absolute difference between current signal and previous stuck signal
DIFF = abs(MV_Percentage_present - xss);
if DIFF > S
% adjust output from valve if DFF is greater than
% (dead-band and stick band), i.e. valve is at
% beginning of trajectory
MV_output = MV_Percentage_present - sign(vnew)*( (S - J)/2 );
else
% keep output at stuck value
MV_output = MV_output_previous;
end
elseif I == 1 % if in stuck condition during moving phase
DIFF = abs(MV_Percentage_present - xss);
if DIFF > J % if DIFF is greater than slip jump
I = 0; % remove flag for being stuck during moving phase
% adjust MV_output
MV_output = MV_Percentage_present - sign(vnew)*( (S - J)/2 );
else
% keep input at the same value
MV_output = MV_output_previous;
end
end
elseif sign(vnew) ~= sign(vold) % if direction of valve movement changes
if sign(vnew) == 0 % if valve reaches a stop position during a transient
I = 1; % indicate this stop condition with flag
xss = MV_Percentage_previous;%MV_incoming_previous; % update memory variable for previous stuck position
MV_output = MV_output_previous; % keep MV output signal constant
else
% do the same for a stuck position not reached during the
% valve's moving phase
xss = MV_Percentage_previous;%MV_incoming_previous;
MV_output = MV_output_previous;
end
end
end % end check for final element saturation
end % end stiction function