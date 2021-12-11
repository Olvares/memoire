function [t,yem] = EulerModif(f,y0,inter,varargin)
%% About this script
% This function solves ODEs and systems of ODEs of the form $y'(t)=f(t,y)$
% using the modified Euler method of the form $y'(t)=f(t,y)$. The implicit
% form of this method is as follows:%
% $y_{i+1} = y_{i} + h*f(t_{i+1},y_{i+1}$ where $y_i = y(t_i)$
%
% Author : Olivier M. ADJAGBA
%
% Last update : 13/07/2021
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
%   yem => solution of the ODE or the ODE system
%
%% Variables
%   y => solution matrix
%

%% ====================== RESOLUTION ======================
% Initialisation
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
    y0 = y0(:);
    y=y0;
    t=inter(1) : h : inter(2);

% Solution
    for i=1:N
        fval = f(t(i),y(:,i));
        y1 = y(:,i) + h*fval; % Explicit Euler explicite calculated for $y_{i+1}$
        fval = f(t(i),y(:,i)) + f(t(i+1),y1);
        y(:,i+1) = y(:,i) + h/2 * fval;
    end
    yem=y';
else
    error(["It is necessary to precise the subdivision step h or the ",...
        "number of time interval subdivision N as {'N',valN} or the ",...
        "time step h as {'h',valh}"])
end
end

