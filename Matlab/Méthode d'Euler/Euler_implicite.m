function [ Ys ] = Euler_implicite( f,y,t,N )
% Cette fonction permet de résoudre les équations différentielles de la
% de y'(t)=f(t,y). La forme implicite de cette méthode est la suivante:
% y(t(i+1))=y(t(i))+h*f(t(i+1),y(t(i+1))
%
h=1/N;
for i=1:N
    y1 = y(i) + h*f(t(i),y(i)); % Euler explicite pour calculer y(t(i+1))
    y(i+1) = y(i) + h*f(t(i+1),y1); % Euler implicite pour recalculer y(t(i+1))
end
Ys=y;
end

