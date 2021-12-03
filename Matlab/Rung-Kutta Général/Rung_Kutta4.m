function [ t, yrk4 ] = Rung_Kutta4(f, y0, inter, N )
%% Cette fonction est un solver des systèmes d'équations différentielles
% d'ordre 1 avec conditions initiales
%% Entrée
% f: est une fonction comportant L'expression de chaque ED d'ordre 1 créé 
% comme dans la fonction diff_Fi.m qui est dans ce dossier.
% y0: est un vecteur colonne comporte la solution initiale.
% n: est l'odre de l'équation différentielle
% inter: est l'intervalle d'etude
% N: est le nombre de subdivision de l'intervalle inter
%% Sortie
% t: est le vecteur comportant les éléments de inter subdivisé en Nh
% sous-intevalles
% yrk4: est la solution du système
%% Variables
% k est la matrice comportant les coéfficients de Rung-Kutta
% y est la matrice solution
% h est le pas de subdivision


%% Résolution

% Initialisation
y0 = y0(:);
y=y0;
h=(inter(2)-inter(1))/N;
t=linspace(inter(1),inter(2),N+1);
n = length(y0);
k=zeros(n,4);

% Algorithme de Rung-Kutta d'ordre 4
for i=1:N
    
    k(:,1)=h*f(t(i),y(:,i));
    
    k(:,2)=h*f(t(i)+h/2,y(:,i)+k(:,1)/2);
    
    k(:,3)=h*f(t(i)+h/2,y(:,i)+k(:,2)/2);
    
    k(:,4)=h*f(t(i)+h,y(:,i)+k(:,3));
    
    y(:,i+1)=y(:,i)+(k(:,1)+2*(k(:,2)+k(:,3))+k(:,4))/6;
    
    %t=t+h;
end
yrk4=y;

end

