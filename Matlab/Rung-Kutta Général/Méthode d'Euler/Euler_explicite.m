function [ Ys ] = Euler_explicite( f,y,t,N )
% Cette fonction permet de résoudre les équations différentielles de la
% de y'(t)=f(t,y). La forme explicite de cette méthode est a suivante:
% y(t(i+1))=y(t(i))+h*f(t(i))
%
h=1/N;
for i=1:N
    y(i+1)=y(i)+h*f(t(i),y(i));
end
Ys=y;
end

