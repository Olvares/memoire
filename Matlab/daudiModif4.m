function [dX] = daudiModif4(t, X, coef, drate, varargin)
%% About this function
% This function is Daudi's model modified by introducing susceptible and
% infected compartiments. Its becomes an epidemiology model SIS in the
% maize population
%
% Author : Olivier M. ADJAGBA
%
% Last update : 12/07/2021
%
%
%% Inputs :
%   X => a vector containing the solution at time $t_i$ in the following order
%   $X = [x_{s_n}n, x_{i_n}, x_n, y_n, z_n, w_n]$ where
%       $x_{s_n}$ => $x$ => maize population at time $t_i$
%       $y_i$ => $y$ => caterpillar population at time $t_i$
%       $z_i$ => $z$ => adult population at time $t_i$
%       $w_i$ => $w$ => egg population at time $t_i$
%   coef => a vector containing the following parameters of Daudi's model
%   in the proposed other : coef = [an, l, d, g, rh, e]
%       an => $\alpha$ or $\eta$ => caterpillar attack rate at vegetative
              % stage ($\alpha$) or at reproductive stage ($\eta$)
%       l => $\lambda$ => maize death due to caterpillar attaque
%       d => $\delta$ => caterpillar to adult rate
%       ga => $\gamma$ => egg to caterpillar
%       rh => $\rho$ => fecondity rate
%       e = $e$ => maize's biomass conversion to caterpillar biomass
%   drate => a vector containing the mrtality rate parameters of Daudi's
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
    k = coef(7);
% Mortality rates
    uy = drate(1);
    uz = drate(2);
    uw = drate(3);

% Others
    if ~isempty(varargin)
        ta = varargin{1}; t0 = varargin{2};
        tw = ta(1) + t0;
        ty = ta(2) + t0;
        b = varargin{3};
    end
%% ========================= VARIABLES ======================
xs = X(1);     % susceptible population of maize
xi = X(2);     % initial population of maize
x = X(3);      % initial population of maize
y = X(4);      % initial population of caterpillard
z = X(5);      % initial population of adult
w = X(6);      % initial population of egg

%% ====================== THE MODEL ========================
%m = @(t,x) l*x;
m = @(t,x) l*exp(-l*t)*x;
%m = @(t,x) l*x*(1-x/k/4);
g = @(x) an*x;
% Possible amelioration
%{
    f = @(x) l*x(1-x/k);
    f = @(x) l*exp(-l*x);
    f = @(x) 0.5*exp(-l*x);
%}

dxs = b*xi - g(xs)*y ;
dxi = g(xs)*y - m(t-t0, xi) - b*xi;
dx = - m(t-t0,xi); % This is the real maize population variation
dy = e*g(xs)*y + (t>=tw)*ga*w - (t>=ty)*d*y - uy*y;
dz = (t>=ty)*d*y - uz*z;
dw = rh*z - (t>=tw)*ga*w - uw*w;

dX = (t>=t0)*[dxs; dxi; dx; dy; dz; dw];
end
