clear all
close all
%clc
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
% Initialisation
n=2; N=20;
t=0; y0=[2;1];
%----------------------------Solution itérative----------------------------
[t1, y] = Rung_Kutta4( @diff_Fi,y0,n,[0 2],N );
%----------------------------Solution analytique---------------------------
ya1=2*exp(t1)-t1.*exp(t1);
ya2=exp(t1)-t1.*exp(t1);
%-----------------------------Représentation-------------------------------
plot(t1,ya1,'*r',t1,ya2,'*b',t1,y(1,:),'b',t1,y(2,:),'r')
grid on
legend('y1(t) analytique','y2 analytique','y1 par Rung-Kutta',...
    'y2 par Rung-Kutta','location','SW')
%% ******************************** FIN ***********************************