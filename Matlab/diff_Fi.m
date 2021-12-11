function [ dif_y ] = diff_Fi( t,y )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dif_y(1,1) = y(2);
dif_y(2,1) = 2*y(2)-y(1);

end

