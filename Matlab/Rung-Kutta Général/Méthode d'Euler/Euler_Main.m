clear all
close all
clc
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Programme de 
%
% Méthode
%
% Critère:
%
% Paramétres
%
% Entrées:
%
% Sortie:
%
% Auteur: ADJAGBA Olivier M.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ******************************* DEBUT **********************************
%-------------------------------Variables----------------------------------
% N,[a,b]

%-----------------------------Initialisation-------------------------------
N=1000;
t=linspace(0,1,N+1);
t1=linspace(0,1,11);
y0=2*pi;
ya=zeros(1,11);
%------------------------------Résolution----------------------------------
%--------------------------Solution analytique-----------------------------
for i=1:11 % y=f(t)
    ya(i)=bissection( @f_Y,y0,5*pi/2,t1(i),eps,N );
end
ymax=bissection(@tx,y0,5*pi/2,1,eps,N); % borne supérieure de y
y2=linspace(2*pi,ymax,11); % décopage de l'intervalle de y
t2=tan(log(tan(y2)+1./cos(y2))); % t=f(y)
%---------------------------Solution itérative-----------------------------
yex=Euler_explicite(@dif_y,y0,t,N); % Algorithme d'Euler explicite
yim=Euler_implicite(@dif_y,y0,t,N); % Algorithme d'Euler implicite
yr=Rung_Kutta4(@dif_y,t,y0,N); % Algorithme de Rung-Kutta
%-----------------------------Représentation-------------------------------
plot(t2,y2,'+k',t1,ya,'*r',t,yex,'g',t,yim,'-.b',t,yr,'-.r')
grid on
legend('Courbe analitique par t(X)','Courbe analitique par X(t)'...
    ,'Courbe d''Euler explicite','Courbe d''Euler implicite',...
    'Courbe de Rung-Kutta','Location','NW');

axis ([0 1 min(yex) max(yex)])
%% ******************************** FIN ***********************************
