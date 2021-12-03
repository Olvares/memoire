function [ x ] = bissection( f,x1,x2,t,epsi,N )
% Cette fonction calcule pour chaque valeur de x par la méthode de la
% bissection
%

for i=1:N
    xm=(x1+x2)/2;
    test=abs((x2-x1)/(2*xm));
    if test<epsi
        break
    elseif f(x1,t)*f(xm,t)<0
        x2=xm;
    elseif f(xm,t)*f(x2,t)<0
        x1=xm;
    end 
end
x=xm;
end

