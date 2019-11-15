function [Bj] = B_coefficient(m,k,c,w,j,Bjp1)
%B_COEFFICIENT Summary of this function goes here
%   Detailed explanation goes here
Bj=(k(j)+c(j)*1i*w)./( -m(j+1)*w.^2+(k(j)+k(j+1))+(c(j)+c(j+1))*1i*w -(k(j+1)+c(j+1)*1i*w).*Bjp1);
end

