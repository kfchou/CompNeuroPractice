% does not work for populations of neurons

dynasimPath = 'DynaSim';
addpath(genpath(dynasimPath));


%% Loading in Spikes into an Input Population
InputNeuron={
    'dV/dt=I; V(0)=-70'
    'if(ismember(t,spiketimes+1))(V=reset)'
    'thresh=-55; reset=-70; I=0; spiketimes=[]'
    
    'if (ismember(t,spiketimes))(I=2000)'
    'if (~ismember(t,spiketimes))(I=0)'
   
    'monitor V.spikes(thresh)'
     };


% define spiketimes
spiketimes = (1:9)'*100 * [1 1];
% dynasim specification
s=[];
s.pops(1).name='IC';
s.pops(1).size=2;
s.pops(1).equations=InputNeuron;
s.pops(1).parameters={'spiketimes',spiketimes};

% simulation
data=dsSimulate(s,'time_limits',[0 1000],'solver','rk1','dt',.01);
dsPlot(data);
dsPlot(data,'plot_type','rastergram');

