function [Dynstiff,Imp ,Dynmass,B] = calculum(m,k,c,w)
%CALCULUM Summary of this function goes here
n=length(m);
B = zeros(length(w),n);

for i=n-1:-1:1
    B(1:length(w),i)=B_coefficient(m,k,c,w,i,B(1:length(w),i+1));
end

Dynstiff=( -m(1)*w.^2+k(1) +1i*w*c(1) - B(:,1).*(k(1)+1i*w*c(1)) );
Imp     =( -m(1)*w.^2+k(1) +1i*w*c(1) - B(:,1).*(k(1)+1i*w*c(1)) )./(1i*w);
Dynmass =Imp./(1i*w);

end

