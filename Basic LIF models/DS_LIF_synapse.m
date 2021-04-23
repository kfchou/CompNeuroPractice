% simulates a LIF neuron in response to a synaptic current with known
% presynaptic spiking

addpath(genpath('C:\Users\Kenny\Desktop\GitHub\DynaSim'));
addpath('C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\mechs');

% time vectors
dt = 0.1;                       % [ms]
t_max = 500;					% [ms]
t_vec = 0:dt:t_max;				% [ms]

% define populations
nCells = 1;
s = struct();
s.populations(1).name='LIF';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {...
    'noise',0,...
    'E_leak',-70,...         % Resting membrane potential [mV]
    'g_leak',1/10,...        % [uS]
    'C',1,...                % membrance capacitance [nF], tau_m=R_m*C_m
    'V_thresh',-54,...       % Threshold for spiking [mV]
    'V_reset',-80,...        % Reset value, after a spike [mV]
    'V_spike',40,...         % Size of a spike [mV]
    'Itonic',0,...         % Injected current [nA]
    'ton',0,...              % [ms]
    'toff',1000,...          % [ms]
    't_ref',0,...            % refractory period [ms]
    };

% Synaptic conductance definition
ICin = {...
    'g_postIC = 0.1'
    'E_exc = 0'
    'ICinput=zeros(1,1000)'
    'I(X,t) = g_postIC .* ICinput(ceil((t+dt)/dt)) .* (X - E_exc);'
    'monitor I'
    '@isyn += I(X,t)'
    };

s.mechanisms(1).name='ICin';
s.mechanisms(1).equations=ICin;

% synaptic conductance due to pre-synaptic spikes
s_tau = 10;                     % Time constant for alpha function [ms]
g_s = 0.024;					% Synaptic conductance [uS]
t_pre = [0.05 0.15 0.19 0.3 0.32 0.4 0.41]*1000;	% Times for synaptic release [ms]
g_vec = 0:dt:100;
g = g_vec/s_tau.*exp(1-g_vec/s_tau); % EPSP waveform
g_post = zeros(1,length(t_vec));
g_post(t_pre*10) = 1;
g_post = conv(g_post,g); % synaptic conductance due to pre-synaptic spikes
g_post(5002:end) = [];

s.connections(1).direction='LIF->LIF';
s.connections(end).mechanism_list='ICin';
s.connections(end).parameters={'g_postIC',g_s,'ICinput',g_post};

t_end = 500; %ms;
data=dsSimulate(s,'dt',0.1,'tspan',[0 t_end]);
for i = 1:length(data)
    data(i).LIF_V(data(i).LIF_V_spikes==1)=40; % insert spikes
end
figure;
subplot(3,1,1); plot(data.time/1000,data.LIF_V/1000); grid;
ylabel('cell voltage [V]')
subplot(3,1,2); plot(data.time/1000,data.LIF_LIF_ICin_I/1E9)
ylabel('current [A]'); xlabel('time (s)'); grid;
subplot(3,1,3); plot(data.time,g_post);

