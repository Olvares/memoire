clear all
close all
clc
%% About this script
% This script is built to solve linear oer non-linear ODE with a numerical
% solver
%
% Methods : RK4, Modified Euler
%
% Author : Olivier M. ADJAGBA
%
% Last update : 12/07/2021
%

%% ================== PARAMETERS' VALUES ========================
% Coefficients
    a = 0.000154;   % $\alpha$ caterpillar attack rate at vegetative stage
    n = 0.000154;   % $\eta$ caterpillar attack rate at reproductive stage
    l = 0.015;      % $\lambda$ maize death due to caterpillar attaque
    d = 0.071;      % $\delta$ caterpillar to adult rate
    g = 0.071;      % $\gamma$ egg to caterpillar
    rh = 0.0417;    % $\rho$ fecondity rate
    k = 500;        % $k$ Maximum number of maize plant in the garden
    e = 1.6;        % $e$ biomass conversion to caterpillar
% Mortality rates
    uy = 0.0071;    % $\mu_y$ mortality rate of caterpillar
    uz = 0.115;     % $\mu_z$ motality rate of adult
    uw = 0.04;      % $\mu_w$ motality rate of egg
    drate = [uy, uz, uw];

%% ======================== TIME ============================
t0 = 0;         % $t_0 = 0$ initial time
t1 = 63;        % $t_1 = 63$ the end of vegetaive stage
T = 160;        % $T = 160$ the end of the reproductive stage
h = 0.1;

%% ====================== RESOLUTION ========================
% figure customization
    fig = 0;
    tit = ["susceptible maize", "infected maize", "maize", "caterpillar", "adult", "egg"];
    xlab = 't';
    ylab = ["x_s(t)", "x_i(t)", "x(t)", "y(t)", "z(t)", "w(t)"];
    colors = 'kbrmgc';
    plotStyle = {'-','--','-.',':.', '*.', '-^'};
    
% Others variables
    numsim = 4; % The number of simulation for diferent value of $z_0$
    nbeq = 6; % number of ODE
    s10 = zeros(numsim, nbeq); % solution at $t = t_1$
    
% Solution for different stage
    for stage = 1:2:3 % 1 correspond to the vegetative stage, 2 for the
                    % reproductive one and three for the whole period
        % Variables initial value
            if stage == 1 % Vegatative stage
                xs0 = k * ones(numsim, 1);
                xi0 = zeros(numsim, 1);
                x0 = xs0;
                y0 = xi0;
                z0 = [1; 30; 45; 60];
                w0 = y0;
                I = [t0 t1]; % time interval
                coef = [a, l, d, g, rh, e, k];
                s0 = [xs0, xi0, x0, y0, z0, w0]; % Matrice for the initial solution
            elseif stage == 2 % Reproductive stage
                coef = [n, l, d, g, rh, e, k];
                s0 = s10; % Matrice for the initial solution from the final
                         % final solution of the vegetative stage
                I = [t1 T]; % time interval
            else % The whole period
                I = [t0 T]; % time interval
                s0 = [xs0, xi0, x0, y0, z0, w0]; % Matrice for the initial solution
            end
            % Intermeadiate function
                interFun = @(t,y) daudiModif(t, y, coef, drate); % To maintain
                    % compatibility with our ODE solvers asked functions
        fig = fig + nbeq*(stage-1);
        for i = 1:numsim
            % [t,s] = RungKutta4(interFun, s0(i,:), I, 'h', h); % RK4
            [t,s] = EulerModif(interFun, s0(i,:), I, 'h', h); % Euler
            s10(i,:) = s(end,:); % Retain the  final solution of the
            % vegetative stage which is also the initial value for the
            % reproductive  stage for the different value of $z_0$

            % for the legend
            name = [sprintf('xs%d = %.0f  ',(stage==2),s0(i,1)),...
                sprintf('xi%d = %.0f  ',(stage==2),s0(i,2)),...
                sprintf('x%d = %.0f  ',(stage==2),s0(i,3)),...
                sprintf('y%d = %.0f  ',(stage==2),s0(i,4)),...
                sprintf('z%d = %.0f  ',(stage==2),s0(i,5)),...
                sprintf('w%d = %.0f  ',(stage==2),s0(i,6))];

            % Solution plot
            for j = 1:nbeq
                figure(fig+j)
                hold on
                % subplot(2,2,j
                if i == numsim
                    xlabel(xlab)
                    ylabel(ylab(j))
                    title(tit(j))
                    legend show
                    legend('location','southoutside')
                    legend('boxoff')
                    grid on
                 end
                plot(t,s(:,j),plotStyle{i}, 'color', colors(i), 'DisplayName', name)
            end
        end
    end