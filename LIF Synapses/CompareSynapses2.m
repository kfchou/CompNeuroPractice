% compare synapses
% one network per synapse
dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
addpath(genpath(dynasimPath));

% basic cell definition
LIF={
    'dV/dt=(E_leak-V+noise*randn-@isyn+R*I(t))/tau; V(0)=-70'
    'if( any(t<tspike + tref, 1) )(V=V_reset)'
    'if(V >= V_thresh)(V = V_reset)'
    'tau=10; tref=1; E_leak=-70; V_thresh=-55; V_reset=-75; R=9; noise=0'
    'Itonic = 0;         % injected current, [nA]'
    'ton = 100;            % [ms]'
    'toff = 900;        % [ms]'
    'I(t)=(Itonic+Itonic*square(2*pi*f*t))*(t>ton&t<toff)'
    'f=0.005'
%     'I(t)=2*(Itonic+Itonic*square(2*pi*f*t))'
    'monitor V.spikes(V_thresh)'
     };

% synapse 1: DS definition
iampa={
  'gSYN=0.5; ESYN=0; tauD=2; tauR=0.4; delay=0'
  'netcon=eye(N_pre,N_post)'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'I(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
  '@isyn += I(V_post)'
  'monitor I'
};

% synapse 2:
dayan={
    'gSYN = 1'
    'ESYN = 0'
    'tauR = 0.4;'
    'tauD = 2;' 
    'delay = 0;'
    'netcon=eye(N_pre,N_post)'
    'epsp(t) =  ( exp(-t/tauD) - exp(-t/tauR) ).*(t>0);'
    'tau = tauR*tauD/(tauR+tauD)'
    'p'' = (epsp(t-tspike_pre(1)-delay)-p)/tau'      % mysterious equation
    'p(0) = 0'
    'I(X,t) = gSYN .* p * netcon .*(X - ESYN)'
    '@isyn += I(V_post,t)'
    'monitor I'
};

%mine
mySyn={
    'gSYN = 1'
    'ESYN = 0'
    'tauR = 0.4;'
    'tauD = 2;' 
    'delay = 0;'
    'netcon=eye(N_pre,N_post)'
    'q(x) = 1*(x<0.5)'
    'alpha(t) = 1/tauR*q(t-tspike_pre(1)-delay)'
    'beta = 1/tauD'
    'm'' = alpha(t)-(alpha(t)+beta)*m'
    'm(0) = 0'
    'I(X,t) = gSYN .* m * netcon .*(X - ESYN)'
    '@isyn += I(V_post,t)'
    'monitor I'
};

% Neuron populations
nPops = 3;
itonic = 1;

ABnetcon = eye(nPops)*1;
ABnetcon(2,2) = 2;
ABnetcon(3,3) = 3;

s=[];
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=iampa;
s.mechanisms(2).name='dayan';
s.mechanisms(2).equations=dayan;
s.mechanisms(3).name='mySyn';
s.mechanisms(3).equations=mySyn;

s.pops(1).name='A';
s.pops(1).size=nPops;
s.pops(1).equations=LIF;
s.pops(1).parameters={'Itonic',itonic};

s.pops(2).name='B';
s.pops(2).size=nPops;
s.pops(2).equations=LIF;

synapses = {s.mechanisms.name};

for i = 1:3
    s.cons(1).direction='A->B';
    s.cons(end).mechanism_list=synapses{i};
    s.cons(end).parameters={'delay',0,'gSYN',3,'netcon',ABnetcon};

    % simulation
    data=dsSimulate(s,'time_limits',[0 1000],'solver','rk1','dt',.01);
    dsPlot(data,'plot_type','rastergram');

    figure;
    subplot(3,1,1);
    plot(data.time, data.A_V); title('A voltage - current pulse input')
    subplot(3,1,2);
    plot(data.time, data.(['B_A_' (synapses{i}) '_I'])); title(['synaptic current:' synapses{i}])
    subplot(3,1,3);
    plot(data.time, data.B_V); title('B voltage')
end
