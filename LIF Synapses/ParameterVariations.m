% 2 neuron network
% one synapse
% compare FR vs parameters
dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
addpath(genpath(dynasimPath));

% parameters
% Neuron populations
nPops = 1;
itonic = 8; % defualt is 1,2, and 3;
gsyn = 2.5; % default is 3;

% vary parameters
varies = struct();
varies(1).conxn = 'A';
varies(end).param = 'Itonic';
varies(end).range = 1:10;

varies(end+1).conxn = 'A->B';
varies(end).param = 'gSYN';
varies(end).range = 0:0.25:7;

varies(end+1).conxn = 'B';
varies(end).param = 'Itonic';
varies(end).range = 0;

varied_param = find(cellfun(@length,{varies.range})>1);
expVar = [varies(varied_param).conxn ' ' varies(varied_param).param];

vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end

% netcons
ABnetcon = eye(nPops)*1;
% ABnetcon(2,2) = 2;
% ABnetcon(3,3) = 3;

% basic cell definition
LIF={
    'tau=10; tref=.5; E_leak=-70; V_thresh=-55; V_reset=-75; R=10; noise=0'
    'dV/dt=( (E_leak-V) + noise*randn -@isyn + R*I(t) )/tau; V(0)=-70'
    'if( any(t<tspike + tref, 1) )(V=V_reset)'
    'if(V >= V_thresh)(V = V_reset)'
    'g_leak = 1/10; C=1;'
    'Itonic = 0;         % injected current, [nA]'
    'ton = 100;            % [ms]'
    'toff = 400;        % [ms]'
    'I(t)=(Itonic+Itonic*square(2*pi*f*t))*(t>ton&t<toff)'
    'f=0.005'
    'monitor V.spikes(V_thresh,2)'
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

s=[];
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=iampa;

s.pops(1).name='A';
s.pops(1).size=nPops;
s.pops(1).equations=LIF;
s.pops(1).parameters={'Itonic',itonic};

s.pops(2).name='B';
s.pops(2).size=nPops;
s.pops(2).equations=LIF;

s.cons(1).direction='A->B';
s.cons(end).mechanism_list='iampa';
s.cons(end).parameters={'delay',0,'gSYN',gsyn,'netcon',ABnetcon,'ESYN',0};

% simulation
data=dsSimulate(s,'vary',vary,'time_limits',[0 400],'solver','rk1','dt',.01);
% dsPlot(data,'plot_type','rastergram');

if isempty(varied_param) == 0
    yspace = 0.95/3;
    height = yspace*0.75;
    width = 0.95/2*0.9;
    xpadSingle = (1-1*width)/2;
    xpadDouble = (1-2*width)/2;
    xstart = 180;
    xend = 320;
    for i = 8
    figure;
    subplot('position',[xpadSingle 0.05+2*yspace width height]);
    plot(data(i).time, data(i).A_V(:,1)); title('A voltage - current pulse input')
    xlim([xstart xend]);
    subplot('position',[xpadSingle 0.05+1*yspace width height]); 
    plot(data(i).time, data(i).B_A_iampa_Isyn(:,1)); title('synaptic current: iampa')
    xlim([xstart xend]);
    subplot('position',[xpadSingle 0.05+0*yspace width height]); 
    plot(data(i).time, data(i).B_V(:,1)); title('B voltage')
    xlim([xstart xend]);
    end
end

% 2d plots
if length(varied_param) == 2
    clear AFR BFR
    nvar1 = length(varies(1).range);
    nvar2 = length(varies(2).range);
    for i = 1:length(data)
        AFR(i) = sum(data(i).A_V_spikes)/0.1;
        BFR(i) = sum(data(i).B_V_spikes)/0.1;
    end
    AFR = reshape(AFR,nvar2,nvar1);
    BFR = reshape(BFR,nvar2,nvar1);
    figure;
    imagesc(varies(1).range,varies(2).range,AFR)
    xlabel('A I tonic')
    ylabel('A->B gSYN')
    title('neuron A FR')

    figure;
    imagesc(varies(1).range,varies(2).range,BFR)
    xlabel('A I tonic')
    ylabel('A->B gSYN')
    title('neuron B FR')

    figure;
    cmap = brewermap(50,'RdBu');
    imagesc(AFR(1,:),varies(2).range,AFR-BFR)
    xlabel('A Firing Rate')
    ylabel('A->B tauD')
    title('Difference in FR (A-B)')
    caxis([-300 300])
    colormap(cmap)
end
% figure;
% plot(varies(varied_param).range,AFR,varies(varied_param).range,BFR)
% xlabel(expVar)
% ylabel('FR (Hz)')
% hl = legend('neuron A','neuron B');
% hl.Location = 'northwest';