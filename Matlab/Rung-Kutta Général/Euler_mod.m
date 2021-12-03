function [t,Ys] = Euler_mod(f,y0,inter,N)
% Cette fonction permet de résoudre les équations différentielles de la
% de y'(t)=f(t,y). La forme implicite de cette méthode est la suivante:
% $y(t(i+1)) = y(t(i)) + h*f(t(i+1),y(t(i+1))$
%
h=1/N;
y0 = y0(:);
y = y0;
h=(inter(2)-inter(1))/N;
t=linspace(inter(1),inter(2),N+1);
for i=1:N
    y1 = y(:,i) + h*f(t(i),y(:,i)); % Euler explicite pour calculer y(t(i+1))
    y(:,i+1) = y(:,i) + h/2 * (f(t(i),y(:,i)) + f(t(i+1),y1));
end
Ys=y;
end

