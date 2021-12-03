function [ yrk ] = Rung_Kutta4( f,t,y,N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

h=1/N;
for i=1:N
    k1=h*f(t(i),y(i));
    k2=h*f(t(i)+h/2,y(i)+k1/2);
    k3=h*f(t(i)+h/2,y(i)+k2/2);
    k4=h*f(t(i)+h,y(i)+k3);
    y(i+1)=y(i)+(k1+2*(k2+k3)+k4)/6;
end
yrk=y;

end

