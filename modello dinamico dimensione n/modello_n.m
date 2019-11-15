%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Modello matematico per i campioni di piccole dimensioni
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

clear all
%close all
load f %vettore f
w=2*pi*f; %vettore omega
n=4; %ordine del sistema
m=ones(n,1);
k=10^9*ones(n,1);
z=ones(n,1);
c=z.*k;
 

[K,I,M]= calculum(m,k,c,w,n,f);

fig3=figure();
plot(f,20*log10(abs(K)))
grid on
fig3.Children.XScale='log';

return

%introduco le stesse costanti del modello di rosario
m1= 1.4293; %massa della piastra di carico pesante 1 [kg];
%m1= 2.8871; %massa della piastra di carico pesante 2 [kg];
%m2= 0.3926; %massa campione 1 [kg]; 
%m2=0.1461; %massa campione 2 [kg];
m2=0.0383; %massa polipropilene
%m3=40; %massa base cemento [kg];
%Costanti elastiche [N/m]: 
kp=59.25e+9;%cost.elast. piastra carico grande [N/m]
%kc2=2.5e7;%cost.elast.AC ref=4e6N/m
E=2*10^9;
s=0.098*0.096;
kc=E*s/0.005
ks=21.6e+9;%cost.elast. base cemento [N/m];
%c=coefficiente di smorzamento=damping ratio*radice quadrata di k su m [Ns/m]: 
%smorzato:c=1,non smorzato:c=0,%fortem smorzato:c>1.
%ref:damping ratio AC=2%-9%;
c1= 0%.00001*2*sqrt(k1*m1); 
c2= 0%.01*2*sqrt(k2*m2);
c3= 0%.03*2*sqrt(k3*m3);    
% m1=1.44;
% m2=30;
m3=15;
k1=kp;
k2=kc;
k3=ks;

w=2*pi*f;
%<<<<<<<<<<<<<<<<<<<<<<<<
% Stampa con masse originali
%<<<<<<<<<<<<<<<<<<<<<<<<
A = (k2+c2*1i.*w)./(-m3 *w.^2 + (k3 + k2) + (c2 + c3)*1i.*w);
B = (k1+c1*1i.*w)./(-m2 *w.^2 + (k1 + k2) + (c1 + c2) *1i.*w - (k2  + c2 *1i.*w).*A);

DMass=(-m1.*w.^2+k2+c1.*1i.*w-(k1+c1.*1i.*w)).*B./(-w.^2);
DSstiff=(-m1.*w.^2+k2+c1.*1i.*w-(k1+c1.*1i.*w)).*B;
I=((-m1.*w.^2+k2+c1.*1i.*w-(k1+c1.*1i.*w)).*B./(1i.*w));
figure (1)
,hold on
p1=plot(f,20*log10(abs(I)), 'r-','LineWidth', 3);


M1=(m1+m2+m3).*w.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M1)), 'r-', 'LineWidth', 1),
set(gca, 'XScale', 'log'),

M2=(m1+m2).*w.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M2)), 'r-', 'LineWidth', 1),
set(gca, 'XScale', 'log'),

M3=m1.*w.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M3)), 'r-', 'LineWidth', 1),
set(gca, 'XScale', 'log'),
set(gca, 'YScale', 'lin'),

k0=115,54;
fa = f(find(20*log10(abs(M1)) > k0,1))
k=4*pi^2*fa^2*m3


% %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% % Modifica masse e definizioni masse efficaci
% %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% 
m3= m3 /2+ m2/2;
m2= m2 /2+ m1 /2;
m1= m1/2;

A = (k2+c2*1i.*w)./(-m3 *w.^2 + (k3 + k2) + (c2 + c3)*1i.*w);
B = (k1+c1*1i.*w)./(-m2 *w.^2 + (k1 + k2) + (c1 + c2) *1i.*w - (k2  + c2 *1i.*w).*A);

DMass=(-m1.*w.^2+k2+c1.*1i.*w-(k1+c1.*1i.*w)).*B./(-w.^2);
DSstiff=(-m1.*w.^2+k2+c1.*1i.*w-(k1+c1.*1i.*w)).*B;
I=((-m1.*w.^2+k2+c1.*1i.*w-(k1+c1.*1i.*w)).*B./(1i.*w));
figure (1)
,hold on
p1=plot(f,20*log10(abs(I)), 'b-.','LineWidth', 3);
%xlim([20, 1e4])



M1=(m1+m2+m3).*w.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M1)), 'b-.', 'LineWidth', 1),

M2=(m1+m2).*w.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M2)), 'b-.', 'LineWidth', 1),

M3=m1.*w.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M3)), 'b-.', 'LineWidth', 1),


Mx= 20*log10(M1)+10*log10(M2./M1);
plot (f, (Mx), 'g-.', 'LineWidth', 1),
k0=115,54;

fa = f(find(20*log10(abs(M1)) > k0,1))
k=4*pi^2*fa^2*m3

set(gca, 'XScale', 'log'),
set(gca, 'YScale', 'lin'),

%%
puppa