function [t, yrk4] = RungKutta4(f, y0, inter, varargin)
%% About this function
% This function is the implementation of RK4 algorithm to resolve ODE
% linear or non-linear system with initial condition
%
% Author : Olivier M. ADJAGBA
%
% Last update : 12/7/2021
%
%% inputs
%   f => a function containing the expression of the equation or the system
%   of equation of order 1 ($dY = f(t,Y)$) to be solved taking necessarily
%   nd only two inputs: the time $t$ and the variable ($Y$) which is a
%   vector in the case of an ODE system and returning only $dY$ which is
%   the value (or the set of values in the case of an ODE system) of the
%   derivative(s) evaluated in $Y$. It returns a column vector
%   y0 => a vector containing the initial solution.
%   inter => the resolution interval
%   varargin => contains two optionals parameters. Those parameters are 
%   N and h where :
%       N => the number of time interval subdivision
%       h => time step
%       It is not necessary to precise N and h together. Only and necessary
%       one of them is needed.
%       The values are intered like this : options = {'N', valN, 'h', valh}
%       with valx the value of the parameter 'x' not  necessary in this 
%       other, but to give a value to one of them.
%       If no value is precised, h=0.1
%
%% Outputs
%   t => a vector containing the elements of inter subdivided into $N$.
%   sub-intervals
%   yrk4 => solution of the ODE or the ODE system
%
%% Variables
%   k => a matrix with Rung-Kutta coefficients
%   y => solution matrix
%   n => the number of ODE
%

%% ====================== RESOLUTION ======================
% Initialisation
    y0 = y0(:);
    y=y0;
    
    if ~isempty(varargin)
        switch varargin{1}
        case 'h'
            h = varargin{2};
        case 'N'
            N = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
        end
    end
    if exist('N','var')
        h=(inter(2)-inter(1))/N;
    elseif exist('h','var')
        N=(inter(2)-inter(1))/h;
    else
        h = 0.1;
        N = (inter(2)-inter(1))/h;
    end

if exist('N','var') && exist('h','var')
    t=linspace(inter(1),inter(2),N+1);
    n = length(y0);
    k=zeros(n,4);

% Solution
    for i=1:N

        k(:,1)=h*f(t(i),y(:,i));

        k(:,2)=h*f(t(i)+h/2,y(:,i)+k(:,1)/2);

        k(:,3)=h*f(t(i)+h/2,y(:,i)+k(:,2)/2);

        k(:,4)=h*f(t(i)+h,y(:,i)+k(:,3));

        y(:,i+1)=y(:,i)+(k(:,1)+2*(k(:,2)+k(:,3))+k(:,4))/6;

        % t=t+h;
    end

    yrk4=y';
    t = t';
else
    error(["It is necessary to precise the subdivision step h or the ",...
        "number of time interval subdivision N as {'N',valN} or the ",...
        "time step h as {'h',valh}"])
end
end
