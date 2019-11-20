function [A] = A_coefficient(m,k,c,w,n)
%A_COEFFICIENT Summary of this function goes here
%   Detailed explanation goes here
A = (k(n-1)+c(n-1)*1i*w)./(-m(n)*w.^2 + (k(n) + k(n-1)) + (c(n-1) + c(n))*1i*w);
end

