function [dX] = daudi(X, coef, drate)
% This function is Daudi's model
%
%% Inputs :
%   X => a vector containing the solution at time $t_i$ in the following order
%   $X = [x_i, y_i, z_i, w_i]$ where
%       $x_i$ => $x$ => maize population at time $t_i$
%       $y_i$ => $y$ => caterpillar population at time $t_i$
%       $z_i$ => $z$ => adult population at time $t_i$
%       $w_i$ => $w$ => egg population at time $t_i$
%   coef => a vector containing the following parameters of Daudi's model
%   in the proposed other : coef = [an, l, d, g, rh, e, k]
%       an => $\alpha$ or $\eta$ => caterpillar attack rate at vegetative
              % stage ($\alpha$) or at reproductive stage ($\eta$)
%       l => $\lambda$ => maize death due to caterpillar attaque
%       d => $\delta$ => caterpillar to adult rate
%       ga => $\gamma$ => egg to caterpillar
%       rh => $\rho$ => fecondity rate
%       k => $k$ => Maximum number of maize plant in the garden at $t=t_0$
%       e = $e$ => maize's biomass conversion to caterpillar biomass
%   drate => a vector containing the mortality rate parameters of Daudi's
%   model in the proposed other : drate = [uy, uz, uw]
%       uy => $\mu_y$ => mortality rate of caterpillars
%       uz => $\mu_z$ => mortality rate of adult moth
%       uw => $\mu_w$ => mortality rate of eggs
%
% Output :
%   dX => a column vector containing values of the ODE evaluated at time
%   $t_{i+1}$
%

%% ================== PARAMETERS' VALUES ========================
% Coefficients
    an = coef(1);
    l = coef(2);
    d = coef(3);
    ga = coef(4);
    rh = coef(5);
    e = coef(6);
% Mortality rates
    uy = drate(1);
    uz = drate(2);
    uw = drate(3);

%% ========================= VARIABLES ======================
x = X(1);      % initial population of maize
y = X(2);      % initial population of caterpillard
z = X(3);      % initial population of adult
w = X(4);      % initial population of egg

%% ====================== THE MODEL ========================
f = @(x) l*x; g = @(x) an*x;

dx = - f(x) - g(x)*y;
dy = e*g(x)*y + ga*w - d*y -uy*y;
dz = d*y - uz*z;
dw = rh*z - ga*w - uw*w;

dX = [dx; dy; dz; dw];
end