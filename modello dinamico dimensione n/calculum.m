function [Dynstiff,Imp ,Dynmass] = calculum(m,k,c,w,n,f)
%CALCULUM Summary of this function goes here

B = zeros(length(f),n-1);
B(:,n-1) = A_coefficient(m,k,c,w,n);

for i=n-2:1
    B=B_coefficient(m,k,c,w,i,B);
end

Dynstiff=( -m(1)*w.^2+k(1) +1i*w*c(1) -B(:,1).*(k(1)+1i*w*c(1)) );
Imp     =( -m(1)*w.^2+k(1) +1i*w*c(1) -B(:,1).*(k(1)+1i*w*c(1)) )./(1i*w);
Dynmass =Imp./(1i*w);

end

