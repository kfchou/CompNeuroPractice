% simulates 2 coupled LIF neurons with either excitatory or inhibitory
% synapses. Both neurons subjected to constant current injections.

dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
mechPath = 'C:\Users\Kenny\Desktop\GitHub\CompNeuroPractice\DynasimMechs';
addpath(genpath(dynasimPath))
addpath(mechPath)

% excitatory vs inhibitory synapse?
excitatory = 1;
if excitatory, ESYN = 0; else ESYN = -80; end

% time vectors
dt = 0.1;                       % [ms]
t_max = 1000;					% [ms]
t_vec = 0:dt:t_max;				% [ms]

% define populations
nCells = 1;
itonic = 4; %1.8 by default
s = struct();
s.populations(1).name='A';
s.populations(end).equations = 'LIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {...
    'noise',0,...
    'E_leak',-70,...         % Resting membrane potential [mV]
    'g_leak',1/10,...        % [uS]
    'C',1,...                % membrance capacitance [nF], tau_m=R_m*C_m
    'V_thresh',-54,...       % Threshold for spiking [mV]
    'V_reset',-80,...        % Reset value, after a spike [mV]
    'V_spike',40,...         % Size of a spike [mV]
    'Itonic',itonic,...      % Injected current [nA]
    'ton',0,...              % [ms]
    'toff',1000,...          % [ms]
    't_ref',0,...            % refractory period [ms]
    'V_init',-70+randn*2  % Initial state
    };

s.populations(2).name='B';
s.populations(end).equations = 'LIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {...
    'noise',0,...
    'E_leak',-70,...         % Resting membrane potential [mV]
    'g_leak',1/10,...        % [uS]
    'C',1,...                % membrance capacitance [nF], tau_m=R_m*C_m
    'V_thresh',-54,...       % Threshold for spiking [mV]
    'V_reset',-80,...        % Reset value, after a spike [mV]
    'V_spike',40,...         % Size of a spike [mV]
    'Itonic',itonic,...      % Injected current [nA]
    'ton',0,...              % [ms]
    'toff',1000,...          % [ms]
    't_ref',0,...            % refractory period [ms]
    'V_init',-70+randn*2  % Initial state
    };

% Synaptic conductance definition
% synAlpha={
%     'gSYN = 0.005; ESYN = 0; tauD = 10; tauR = 0.5; delay = 0'
%     'netcon=ones(N_pre,N_post)'
%     'epsp(t) = (exp(-t/tauD)).*(t>0)'
%     'p'' = (exp(1)*epsp(t-tspike_pre-delay)-p)/tauD'
%     'p(0) = 0'
%     'I(X,t) = gSYN .* p * netcon .*(X - ESYN)'
%     '@isyn += I(V_post,t)'
%     'monitor functions'
%   };

% dynasim examples (iAMPA, etc) use the following synapse; it does not
% yield the same result.
synAlpha={
    'gSYN = 1; ESYN = 0; tauD = 10; tauR = 0.5; delay = 0'
    'netcon=ones(N_pre,N_post)'
    'epsp(t) = (exp(-t/tauD)).*(t>0)'
    'I(X,t) = gSYN .* sum(epsp(t-tspike_pre-delay)) * netcon .*(X - ESYN)'
    '@isyn += I(V_post,t)'
    'monitor functions'
  };

s.mechanisms(1).name = 'synAlpha';
s.mechanisms(1).equations = synAlpha;

s.connections(1).direction='A->B';
s.connections(end).mechanism_list='synAlpha';
s.connections(end).parameters={'tauD',10,...        % [ms]
                               'gSYN',0.005,...     % [uS]
                               'ESYN',ESYN,...
                              };
                          
s.connections(2).direction='B->A';
s.connections(end).mechanism_list='synAlpha';
s.connections(end).parameters={'tauD',10,...        % [ms]
                               'gSYN',0.005,...     % [uS]
                               'ESYN',ESYN,...
                              };                          

data=dsSimulate(s,'dt',0.1,'tspan',[0 t_max]);
for i = 1:length(data)
    data(i).A_V(data(i).A_V_spikes==1)=40; % insert spikes
    data(i).B_V(data(i).B_V_spikes==1)=40; % insert spikes
end
figure;
subplot(4,1,1); plot(data.time/1000,data.A_V/1000);
ylabel('cell 1 voltage [V]')
subplot(4,1,2); plot(data.time/1000,data.B_V/1000);
ylabel('cell 2 voltage [V]')
subplot(4,1,3); plot(data.time/1000,data.A_B_synAlpha_I);
ylabel('current [A]'); xlabel('time (s)');
subplot(4,1,4); plot(data.time/1000,data.B_A_synAlpha_I);
ylabel('current [A]'); xlabel('time (s)');

figure;
subplot(4,1,1); plot(data.time/1000,data.A_V/1000);
ylabel('cell 1 voltage [V]')
subplot(4,1,2); plot(data.time/1000,data.B_V/1000);
ylabel('cell 2 voltage [V]')
subplot(4,1,3); plot(data.time/1000,data.A_B_synAlpha_I);
ylabel('current [A]'); xlabel('time (s)');
subplot(4,1,4); plot(data.time/1000,data.B_A_synAlpha_I);
ylabel('current [A]'); xlabel('time (s)');
for i = 1:4
    subplot(4,1,i);
    xlim([0.9 1.0]);
end

figure; 
subplot(2,1,1); plot(data.time,data.A_B_synAlpha_p)
ylabel('conductance trace, s_1')
subplot(2,1,2); plot(data.time,data.A_B_synAlpha_epsp)
ylabel('conductance trace, s2_1')
