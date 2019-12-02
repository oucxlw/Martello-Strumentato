freq=0.1:0.01:10^4;

m=50; %m=1.4*10^8;
k=10^2; %k=2.8*10^7;
% Zeta=1./(2*100*pi*freq);
% c=Zeta*k;
c=0.01;

[radicale,phi,angolo,w] = radice(m,k,c,freq);

% fig=figure();
% plot(w,phi)
% hold on
% plot(w,angolo,'--')
% fig.Children.XScale='log';
%angolo è dewrappato rispetto a phi

Impedenza=(k*radicale./w).*exp(1i*(phi+pi/2));

fig2=figure();
plot(freq,20*log10(abs(Impedenza)))
fig2.Children.XScale='log';
grid on
hold on

% valori sospettati essere la causa della risonanza a bassissima frequenza
% m3=50;
% k3=10^2;
% c3=0.001;

m=[1.4 0.8 50];
k=[5*10^10 2.8*10^7 10^2]; %k=[5*10^10 2.8*10^8 10^2];
c=[0.001 0.001 0.001];

A=( k(2)+1i*w*c(2) )./(-m(3)*w.^2+k(2)+k(3)+1i*w.*(c(2)+c(3)));
B=( k(1)+1i*w*c(1) )./(-m(2)*w.^2+k(1)+k(2)+1i*w.*(c(1)+c(2))-A.*(k(2)+1i*w.*c(2)));
Dynstiff=( -m(1)*w.^2+k(1) +1i*w*c(1) -B.*(k(1)+1i*w*c(1)) );
Imp=( -m(1)*w.^2+k(1) +1i*w*c(1) -B.*(k(1)+1i*w*c(1)) )./(1i*w);

fig3=figure();
plot(freq,20*log10(abs(Imp)))
grid on
fig3.Children.XScale='log';

Imp_resample=interp1(freq,Imp,f);