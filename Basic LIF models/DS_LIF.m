dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
mechPath = 'C:\Users\Kenny\Desktop\GitHub\CompNeuroPractice\DynasimMechs';
addpath(genpath(dynasimPath))
addpath(mechPath)

% simulates a LIF neuron in response to a constant current injection

nCells = 1;

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
    'Itonic',1.8,...         % Injected current [nA]
    'ton',0,...              % [ms]
    'toff',1000,...          % [ms]
    't_ref',0,...            % refractory period [ms]
    };

I_e = linspace(0,4,16); % [nA]
vary={
  'LIF','Itonic',I_e;     % amplitude of tonic input to E-cells
  };

t_end = 1000; %ms;
data=dsSimulate(s,'dt',0.1,'tspan',[0 t_end],'vary',vary);
for i = 1:length(data)
    data(i).LIF_V(data(i).LIF_V_spikes==1)=40; % insert spikes
    spkCount(i) = sum(data(i).LIF_V_spikes);
end

% calculate & plot FI curve
figure;
p = cellpairs2struct(s.populations(1).parameters);
taum = p.C/p.g_leak*1E-3; % convert to seconds
fr_analyt = 1./( taum *log(p.E_leak+1/p.g_leak*I_e-p.V_reset) - taum *log(p.E_leak+1/p.g_leak*I_e-p.V_thresh) );
fr_analyt(find(imag(fr_analyt)~=0)) = 0;
plot(I_e*1E-9,spkCount,'o',I_e*1E-9,fr_analyt,'x'); hold on;
xlabel('Injected current [nA]');
ylabel('Firing Rate [Hz]');

relu = max((I_e-1),0)*45;
plot(I_e*1E-9,relu)

legend('numerical','analytica','relu approx')