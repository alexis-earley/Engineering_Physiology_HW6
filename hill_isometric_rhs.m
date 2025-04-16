%hill_rhs.m
% right-hand side for isotonic Hill model
% input x = current length, in cm (scalar)
% output = dx/dt
function dFdt = hill_isometric_rhs(t,F)
global kpe kse b xstar x delay A tv

% scaling function for A(t) as a function of length
s = (x>0.5*xstar)*(x<1.5*xstar)*(cos(pi*(x-xstar)));

% instead of calculating A(t), use precalculated version
%Aloc = 1*(t>delay)*(48144*exp(-(t-delay)/0.0326) - 45845*exp(-(t-delay)/0.034));
% imin = index of nearest match to current value of t
% using a fine grid, no need to interpolate
[tmin,imin]=min(abs(tv-t));
Aloc = A(imin);

dFdt = kse/b*(kpe*(x-xstar)*(x>xstar) - (1+kpe/kse)*F + Aloc*s);

