%% 1. original dynasim example
LIF={
    'dV/dt=(E-V+R*I+noise*randn-@isyn)/tau; V(0)=-65'
    'if( any(t<tspike + tref, 1) )(V=reset)'
    'tau=10; tref=10; E=-70; thresh=-55; reset=-75; R=9; I=1.55; noise=100'
    'monitor V.spikes(thresh,2)'
     };

iampa={
  'gSYN=.5; ESYN=0; tauD=2; tauR=0.4; delay=15'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(X) = gSYN.*sum(f(t-tspike_pre-delay)).*(X-ESYN)'
  '@isyn += Isyn(V_post)'
  'monitor Isyn'
};

s=[];
s.populations(1).name='E';
s.populations(1).equations=LIF;
s.populations(2).name='I';
s.populations(2).equations=LIF;
s.populations(2).parameters={'I',0,'noise',0};
s.connections(1).direction='E->I';
s.connections(1).mechanism_list='iampa';
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=iampa;

data=dsSimulate(s,'time_limits',[0 200], 'solver','rk1', 'dt',.01);
dsPlot(data);
figure;
subplot(4,1,4); plot(data.time,ones(1,length(data.time))*1.55); ylabel({'input current (5khz)','half-wave rectified'}); xlabel('time (ms)')
subplot(4,1,3); plot(data.time,data.E_V); ylabel('presynaptic voltage')
subplot(4,1,2); plot(data.time,data.I_E_iampa_Isyn); ylabel('synaptic current')
subplot(4,1,1); plot(data.time,data.I_V); ylabel('postsynaptic voltage')

%% 2. input current as sinusoid
LIF={
    'dV/dt=(E-V+R*Ie(ceil((t+dt)/dt))+noise*randn-@isyn)/tau; V(0)=-65'
    'if( any(t<tspike + tref, 1) )(V=reset)'
    'tau=10; tref=1; E=-70; thresh=-55; reset=-75; R=9; Ie=zeros(1,1000); noise=100'
    'monitor V.spikes(thresh,2)'
     };

iampa={
  'gSYN=.5; ESYN=0; tauD=2; tauR=0.4; delay=0'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(X) = gSYN.*sum(f(t-tspike_pre-delay)).*(X-ESYN)'
  '@isyn += Isyn(V_post)'
  'monitor Isyn'
};

dt = 0.01; %ms
t_end = 50; %ms
t = 0:dt:t_end; %ms
freq = 5; %khz
% input = cos(2*pi*freq*t)*1000;
% input(input<0) = 0;
input = ones(1,length(t))*20;

s=[];
s.populations(1).name='E';
s.populations(1).equations=LIF;
s.populations(1).parameters={'Ie',input,'noise',0};

s.populations(2).name='I';
s.populations(2).equations=LIF;
s.populations(2).parameters={'Ie',zeros(size(input)),'noise',0};

s.populations(3).name='A';
s.populations(end).equations=LIF;
s.populations(end).parameters={'Ie',zeros(size(input)),'noise',0};

s.connections(1).direction='E->I';
s.connections(1).mechanism_list='iampa';
s.connections(1).parameters={'gSYN',4};

s.connections(end+1).direction='I->A';
s.connections(end).mechanism_list='iampa';
s.connections(end).parameters={'gSYN',4};

s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=iampa;

data=dsSimulate(s,'time_limits',[0 t_end], 'solver','rk1', 'dt',dt);
figure;
subplot(4,1,4); plot(data.time,input); ylabel({'input current (5khz)','half-wave rectified'}); xlabel('time (ms)')
subplot(4,1,3); plot(data.time,data.E_V); ylabel('presynaptic voltage')
subplot(4,1,2); plot(data.time,data.I_E_iampa_Isyn); ylabel('synaptic current')
subplot(4,1,1); plot(data.time,data.I_V); ylabel('postsynaptic voltage')
