function [t,s] = myMod(s0, coef, drate, I)
% This function is built to solve my model with a numerical ODE solver
% denoted odeSolver
%

%% Parameters
an = coef(1);   % $\alpha$n caterpillar attaque rate
l = coef(2);    % $\lambda$ maize death due to caterpillar attaque
d = coef(3);    % $\delta$ caterpillar to adult rate
ga = coef(4);   % $\gamma$ egg to caterpillar
rh = coef(5);   % $\rho$ fecondity rate
e = coef(6);    % $e$ caterpillar grow due to maize consumption
k = coef(7);    % $k$ maize population limite
rr = coef(8);   % $rr$ maize resistace rate
uy = drate(1);  % $\mu_y$ mortality rate of caterpillar
uz = drate(2);  % $\mu_z$ motality rate of adult
uw = drate(3);  % $\mu_w$ motality rate of egg

%% Variable
x0 = s0(1);      % initial population of maize
y0 = s0(2);      % initial population of caterpillard
z0 = s0(3);      % initial population of adult
w0 = s0(4);      % initial population of egg

%% Résolution
f = @(x) - l*x;%*(x/k - 1); %f(y(1))
g = @(x) an*x ;
myf = @(t,y) [f(y(1)) - g(y(1))*y(2);% * exp(-rr*t);
    - e*g(y(1))*y(2) + ga*y(4) - d*y(2) - uy*y(2);
    e*g(y(1))*y(2) + d*y(2) - uz*y(3);
    rh*y(3) - ga*y(4) - uw*y(4)];
N = 50;
[t,s] = Rung_Kutta4(myf, s0, I, I(end)*N);
end