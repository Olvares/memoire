function [t,s] = myModSIS_6(s0, coef, drate, ta, t0, sig, I)
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
%sig = coef(9);  % $\sigma$ immigratin rate

uy = drate(1);  % $\mu_y$ mortality rate of caterpillar
uz = drate(2);  % $\mu_z$ motality rate of adult
uw = drate(3);  % $\mu_w$ motality rate of egg

twy = ta(1)+t0 ;
tyz = twy + ta(2) ;

%% Variable
%x0 = s0(1);      % initial population of maize
%y0 = s0(2);      % initial population of caterpillard
%z0 = s0(3);      % initial population of adult
%w0 = s0(4);      % initial population of egg

%% Résolution
g = @(x) an*x;

m = @(t,xi) l*xi; % Mortality due to cla
%m = @(t,xi) exp(-l * (t-t0))*xi;%exppdf(t, 1/l); % Mortality due to cla
%m = @(t,xi) l*exp(-l * (t-t0))*xi;%exppdf(t, 1/l); % Mortality due to cla
A = @(t,xs,x,y) g(xs) * y / x;%(xs+xi); % Attack
B = @(xi) rr * xi; % Retablishment
Dyz = @(t,xi,x,y) ((t-tyz)>=0) * d * xi * y / x;%(xs+xi); % Caterpillar to Adult
R = @(z) rh * z; %* (1 - z/k); % Reproduction
Dwy = @(t,w) ((t-twy)>=0) * ga * w; % Egg to Larvae

myfSIS = @(t,y) (t>=t0) * [B(y(2)) - A(t,y(1),y(3),y(4));
    A(t,y(1),y(3),y(4)) - B(y(2)) - m(t,y(2));
    -m(t,y(2));% * y(4);
    Dwy(t,y(6)) - Dyz(t,y(2),y(3),y(4)) - (uy + 0.0*y(4))*y(4);
    Dyz(t,y(2),y(3),y(4)) - uz*y(5) + sig;
    R(y(5)) - Dwy(t,y(6)) - uw*y(6)
];
N = 50;
[t,s] = Rung_Kutta4(myfSIS, s0, I, I(end)*N);
%[t,s] = ode45(myfSIS, I, s0);
end