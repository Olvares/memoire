clear all
close all
clc
%% Daudi's model
%% Variables
    x0 = 500;
    y0 = 0;
    z0 = [15, 30, 45 60];
    w0 = 0;
    ta = 10;
%% Parameters
% Coefficients
    an = 0.000154;  % $\alpha$ caterpillar attack rate
    l = 0.015;%0;%.001;%0.015;      % $\lambda$ maize death due to caterpillar attaque
    d = 0.071;      % $\delta$ caterpillar to adult rate
    g = 0.071;      % $\gamma$ egg to caterpillar
    rh = 0.0417;%0.1,0.0417;    % $\rho$ fecondity rate
    e = 1.6; %0.1;%1.6       % $e$ biomass conversion to caterpillar
    k = x0/2;
    coef = [an, l, d, g, rh, e, k];
% Mortality rates
    uy = 0.0071;    % $\mu_y$ mortality rate of caterpillar
    uz = 0.115;     % $\mu_z$ motality rate of adult
    uw = 0.04;      % $\mu_w$ motality rate of egg
    drate = [uy, uz, uw];
% Time
    t0 = 0;
    t1 = 63;
    T = 160;
    I = [t0 T];
%% Resolution
fig = 0;
tit = ["maize", "caterpillar", "adult", "egg"];
for i = 1:length(z0)
    s0 = [x0, y0, z0(i), w0]; % Initial solution with z0 defent value
    %[t,s] = daudiMod(s0, coef, drate, I); % Solution with RK4 method
    [t,s] = daudiMod(s0, coef, drate, I, "EM"); % Solution with Euler method
    for j = 1:4
        figure(fig+j)
        hold on
        %subplot(2,2,j)
        plot(t,s(j,:))
        title(tit(j))
        grid minor
    end
end
%figure(1)
%hold on
%plot(t, x0 * exp(-l*t)) % for verification
%fig = fig + j+1;