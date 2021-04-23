% compare synapses
% mixing synapses within one network
dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
addpath(genpath(dynasimPath));

% parameters
% Neuron populations
nPops = 1;
itonic = 8; % defualt is 1,2, and 3;
gsyn = 3; % default is 3;

ABnetcon = eye(nPops)*1;
% ABnetcon(2,2) = 2;
% ABnetcon(3,3) = 3;
BCnetcon = ABnetcon;
BDnetcon = ABnetcon;

% basic cell definition
LIF={
    'dV/dt=( (E_leak-V) + noise*randn -@isyn + R*I(t) )/tau; V(0)=-70'
    'if( any(t<tspike + tref, 1) )(V=V_reset)'
    'if(V >= V_thresh)(V = V_reset)'
    'tau=10; tref=.5; E_leak=-70; V_thresh=-55; V_reset=-75; R=10; noise=0'
    'g_leak = 1/10; C=1;'
    'Itonic = 0;         % injected current, [nA]'
    'ton = 100;            % [ms]'
    'toff = 900;        % [ms]'
    'I(t)=(Itonic+Itonic*square(2*pi*f*t))*(t>ton&t<toff)'
    'f=0.005'
    'monitor V.spikes(V_thresh,3)'
     };

% synapse 1: DS definition
iampa={
  'gSYN=0.5; ESYN=0; tauD=2; tauR=0.4; delay=0'
  'netcon=eye(N_pre,N_post)'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
  '@isyn += Isyn(V_post)'
  'monitor Isyn'
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
    'p'' = (sum(epsp(t-tspike_pre-delay))-p)/tau'      % mysterious equation
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
    'alpha(t) = sum(1/tauR*q(t-tspike_pre-delay))'
    'beta = 1/tauD'
    'm'' = alpha(t)-(alpha(t)+beta)*m'
    'm(0) = 0'
    'I(X,t) = gSYN .* m * netcon .*(X - ESYN)'
    '@isyn += I(V_post,t)'
    'monitor I'
};


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

s.pops(end+1).name='C';
s.pops(end).size=nPops;
s.pops(end).equations=LIF;

s.pops(end+1).name='D';
s.pops(end).size=nPops;
s.pops(end).equations=LIF;

s.cons(1).direction='A->B';
s.cons(end).mechanism_list='iampa';
s.cons(end).parameters={'delay',0,'gSYN',gsyn,'netcon',ABnetcon};

s.cons(end+1).direction='B->C';
s.cons(end).mechanism_list='mySyn';
s.cons(end).parameters={'delay',0,'gSYN',gsyn,'netcon',BCnetcon};

s.cons(end+1).direction='B->D';
s.cons(end).mechanism_list='dayan';
s.cons(end).parameters={'delay',0,'gSYN',gsyn,'netcon',BDnetcon};

% simulation
data=dsSimulate(s,'time_limits',[0 1000],'solver','rk1','dt',.01);
dsPlot(data,'plot_type','rastergram');

yspace = 0.95/5;
height = yspace*0.75;
width = 0.95/2*0.9;
xpadSingle = (1-1*width)/2;
xpadDouble = (1-2*width)/2;
xstart = 180;
xend = 320;
for i = 1:nPops
figure;
subplot('position',[xpadSingle 0.05+4*yspace width height]);
plot(data.time, data.A_V(:,i)); title('A voltage - current pulse input')
xlim([xstart xend]);
subplot('position',[xpadSingle 0.05+3*yspace width height]); 
plot(data.time, data.B_A_iampa_Isyn(:,i)); title('synaptic current: iampa')
xlim([xstart xend]);
subplot('position',[xpadSingle 0.05+2*yspace width height]); 
plot(data.time, data.B_V(:,i)); title('B voltage')
xlim([xstart xend]);
subplot('position',[xpadDouble 0.05+1*yspace width height]); 
plot(data.time, data.C_B_mySyn_I(:,i)); title('synaptic current: Analytical')
xlim([xstart xend]);
subplot('position',[xpadDouble 0.05+0*yspace width height]); 
plot(data.time, data.C_V(:,i)); title('C voltage')
xlim([xstart xend]);

subplot('position',[xpadDouble+width*1.1 0.05+1*yspace width height]); 
plot(data.time, data.D_B_dayan_I(:,i)); title('synaptic current: Dayan&Abbott')
xlim([xstart xend]);
subplot('position',[xpadDouble+width*1.1 0.05+0*yspace width height]); 
plot(data.time, data.D_V(:,i)); title('D voltage')
xlim([xstart xend]);
end
