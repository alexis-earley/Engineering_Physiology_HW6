% hill_isometric.m
% Simulation of Hill model
% with thanks to Reza Shademer
% see http://www.bme.jhu.edu/~reza/book/muscle_model/musclemodel.htm
% isometric version -- x held constant
% this model accounts for static nonlinearity of the contractile element
% force A with changing length, but does not account for nonlinearities of
% the series and parallel elastic elements (springs)

clear
%close all
global kpe kse b xstar x delay A tv

kpe = 75;   % spring constant of parallel element, g/cm
kse = 136;  % spring constant of series element, g/cm
b   = 50;   % viscosity of parallel dashpot, (g*s)/cm
delay = 0.1; % delay before stimulation, s

xv    = [0.51 0.75 1 1.25 1.49];     % fixed length, cm
xstar = 1;     % resting length, cm
tspan = [0 1];  % time span, s

% pre-calculate the activation function A(t)
dtv = 0.0001;
tv = [tspan(1):dtv:tspan(2)];
A = 1*(tv>delay).*(48144*exp(-(tv-delay)/0.0326) - 45845*exp(-(tv-delay)/0.034));
% generate a train of activation pulses at frequency f
% f=20; pt=(1/f:1/f:1); 
% for n = 1:length(pt), A = A + 1*(tv>delay+pt(n)).*(48144*exp(-(tv-delay-pt(n))/0.0326) - 45845*exp(-(tv-delay-pt(n))/0.034)); end;
cvect ='bgrck';
figure
hold;
for i=1:length(xv)
    x = xv(i);
    F0 = (x>xstar)*kpe*(x-xstar)/(1+kpe/kse);  %re-calculate initial force = steady-state force
    [t,F]=ode15s('hill_isometric_rhs',tspan,F0);
    plot(t,F,cvect(i))
end

