dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
addpath(genpath(dynasimPath));

%% Interacting E and I populations
% tau [ms]: membrane time constant (RC)
% tabs [ms]: absolute refractory period
% delay [ms]: axonal delay
% tspike: reserved variable for storing past spike times when monitoring spikes

% changes:
% Isyn(V_post) -> Isyn(X)
%
LIF={
    'dV/dt=(E-V+(R*Iapp(t)+noise*randn(1,N_pop))+@isyn)/tau; V(0)=-70'
    'if(any(t<tspike+tabs,1))(V=reset)'
    'tau=10; tabs=1; E=-70; thresh=-55; reset=-75; R=9; I=0; noise=0;'
    
    % applied current waveform
    'mask=ones(1,N_pop); f=0.005'; % input mask and input frequency in kHz
    'Iapp(t)=mask.*(I+I*square(2*pi*f*t))';
    
    'monitor Iapp'
    'monitor V.spikes(thresh)'
     };
 
iampa={
  'gSYN=.5; ESYN=0; tauD=2; tauR=0.4; delay=15'
  'netcon=ones(N_pre,N_post)'
  'f(x) = 1*(exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
  '@isyn += -Isyn(V_post)'
  'monitor Isyn'
};

% input mask
Ne=3; Ni=3;
maskE=zeros(1,Ne);
maskE(1)=1; % scale factor for injected current
maskE(2)=2;
maskE(3)=3;

%connectivity matrix
Npre=Ne; Npost=Ni; Nmax=max(Npost,Npre);
netconEI=eye(Nmax);

% dynasim specification
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=LIF;
s.pops(1).parameters={'mask',maskE,'I',1};
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=LIF;
s.cons(1).direction='E->I';
s.cons(1).mechanism_list='iampa';
s.cons(1).parameters={'netcon',netconEI,'delay',0,'gSYN',3};
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=iampa;

% simulation
data=dsSimulate(s,'time_limits',[0 1000],'solver','rk1','dt',.01);
dsPlot(data);
dsPlot(data,'plot_type','rastergram');

figure;
for i = 1:3
subplot(3,1,i)
plot(data.time, data.I_E_iampa_Isyn(:,i))
end
