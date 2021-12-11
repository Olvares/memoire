function [t,s] = myModSIS_2(s0, coef, drate, I)
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
%x0 = s0(1);      % initial population of maize
%y0 = s0(2);      % initial population of caterpillard
%z0 = s0(3);      % initial population of adult
%w0 = s0(4);      % initial population of egg

%% Résolution
f = @(x) - l*x;%*(x/k - 1); f(y(1))
m = @(t) l*exp(-l*t);%exppdf(t, 1/l);
g = @(x) an*x;
myfSIS = @(t,y) [rr * y(2) - g(y(1))*y(4);%/s0(1);% * exp(-rr*t);
    - rr * y(2) + g(y(1))*y(4) - m(t) * y(2);% * y(4);
    - m(t) * y(2);% * y(4);
    - e * g(y(1))*y(4) + ga*y(6) - d*y(4) - uy*y(4);
    (e*g(y(1)) + d) * y(4) - uz*y(5);
    rh*y(5) - ga*y(6) - uw*y(6)];
N = 50;
[t,s] = RungKutta4(myfSIS, s0, I, 'N', I(end)*N);
end