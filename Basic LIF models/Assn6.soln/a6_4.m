% Assignment #6, Per Jesper Sj�str�m% 2/28/00, The I&F neuron% This program has a synaptic conductance% The synaptic conductance is governed by an alpha function%%% Init% CellV_rest = -70e-3;				% The resting membrane potential [V]R_m = 10e6;						% Input resistance of the cell [Ohm]tau_m = 10e-3;					% Membrane time constant [s]V_th = -54e-3;					% Threshold for spiking [V]V_reset = -80e-3;				% Reset value, after a spike [V]V_max = 40e-3;					% Size of a spike [V]V = V_rest;						% Initiate the membrane potential to the leakage reversal [V]% Synapses_tau = 0.01;					% Time constant for alpha function [s]s_max = 1;						% Max amplitude of synaptic gating [unitless]g_s = 0.025e-6;					% Synaptic conductance [S]E_s = 0;						% Synaptic reversal potential [V]s = 0;							% Synapse first order gating is closed [unitless]s2 = 0;							% Synapse second order gating is closed [unitless]t_rel_vec = [0.05 0.15 0.19 0.3 0.32 0.4 0.41];	% Times for synaptic release [s]e = exp(1);						% The number ei_syn = 0;						% Synaptic current [A]% Vectorsdt = 0.0001;					% [s]t_max = 0.5;					% [s]t_vec = 0:dt:t_max;				% [s]V_vec = zeros(1,length(t_vec));	% [V]i_syn_vec = zeros(1,length(t_vec));	% synaptic current [A]%%% Loopfor i = 1:length(t_vec),	V_vec(i) = V;				% Store away the membrane voltage	i_syn_vec(i) = i_syn;		% Store away the synaptic current		if (V==V_max)				% Previous step was a spike?		V = V_reset;			% If so, do a reset	end;	if (length(find(abs(t_rel_vec-t_vec(i))<eps))>0)		s2 = s_max;	end;	s = s + dt*((e*s2-s)/s_tau);	s2 = s2 + dt*(-s2/s_tau);    s1(i) = s;    s21(i) = s2;	tau_v = tau_m/(1+R_m*g_s*s);	V_inf = (V_rest+R_m*g_s*s*E_s)/(1+R_m*g_s*s);	i_syn = -g_s*s*(E_s-V);	V = V_inf + (V-V_inf)*exp(-dt/tau_v);	if (V>V_th),				% Above threshold?		V = V_max;				% If so, produce the spike	end;end;%%% Plotfigure(1);orient tall;subplot(3,1,1);plot(t_vec,V_vec);grid;xlabel('time [s]');ylabel('ampl [V]');subplot(3,1,2);plot(t_vec,i_syn_vec);hold on;for i = 1:length(t_rel_vec),	plot([t_rel_vec(i) t_rel_vec(i)],[0 min(i_syn_vec)],'--r');end;hold off;grid;xlabel('time [s]');ylabel('current [A]');subplot(3,1,3);plot(t_vec,s1); ylabel('pre-synaptic dependent condutance');