% 02) Modello massa-molla-smorzatore
% 3 masse (piastra di carico, provino, base cemento)
% 3 molle (piastra di carico, provino, gommine sotto base cemento)
% 3 smorzatori (piastra di carico, provino, gommine sotto base cemento)
%dovrebbe valere solo per "non attaccato"
% trova Dstiff=F/X1

set (0,'DefaultFigureWindowStyle','docked')
clc
%%
% N.B.: Dopo aver aperto la figura "Collezione Dstiff-....", applico modello
% ad es.: Collezione Dstiff-CCampione1, non attaccato, di lato, rip1-Acc_0-auto-M-cilindrica pesante 1.fig

% Importo i dati di input
clear variables
Inputs=load('Outputs.mat');
f = Inputs.frequenze.f;
f_fft = Inputs.frequenze.f_fft;
sng = Inputs.sng;
win_1 = Inputs.win_1;
PSD = Inputs.PSD;
FFT = Inputs.FFT;
win_A = Inputs.win_A;
result = Inputs.result;
conf = Inputs.conf;

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf,piastre);

m = piastre.massa(conf.piastra); %massa della piastra in uso
h = campioni.h(conf.campione);
s = pi*(campioni.d(conf.campione)/2)^2;

clear Inputs;
omega=2*pi.*f;
m1= piastre.massa(conf.piastra); %massa della piastra di carico pesante 1 [kg];
%m1= 2.8871; %massa della piastra di carico pesante 2 [kg];
m2= campioni.massa(conf.campione);; %massa campione 1 [kg]; 
%m2=0.1461; %massa campione 2 [kg];
m3=15; %massa base cemento [kg] 

m3= m3 /3+2*m2 /3;
m2= m2 /3+2*m1 /3;
m1= m1/3;

%Costanti elastiche [N/m]: 
% k1=2e9;%cost.elast. piastra carico grande [N/m]
% k2=2.5e7;%cost.elast.AC ref=4e6N/m
% k3=4e8;%cost.elast. base cemento [N/m];

k1=40e9;%cost.elast. piastra carico grande [N/m]
k2=2.8e+9;%cost.elast.AC ref=4e6N/m
k3=21.6e+9;%cost.elast. base cemento [N/m];

%c=coefficiente di smorzamento=damping ratio*radice quadrata di k su m [Ns/m]: 
%smorzato:c=1,non smorzato:c=0,%fortem smorzato:c>1.
%ref:damping ratio AC=2%-9%;
c1= 0.00001*2*sqrt(k1*m1); 
c2= 0.01*2*sqrt(k2*m2);
c3= 0.03*2*sqrt(k3*m3);    


h=0.05; %altezza slab AC [m];
S=(pi*0.004)^2; %area impatto martello [m^2]
s=1i.*omega;
s2=-omega.^2;
s3=(1i^3).*omega.^3;
s4=omega.^4;

B=(m2*m3).*s4 + (m2*(c2+c3)+m3*(c1+c3)).*s3 + (c1*(c2+c3)+c2*c3+m2*(k2+k3)+m3*(k1+k2)).*s2...
    +(k2*(c1+c2)+k3*(c1+c2)+k1*(c2+c3)+k2*(c2+c3)-2*c2).*s+(k1*(k2+k3)+k2*k3);
C=(-m3*c1).*s3-(c1*(c2+c3)+m3*k1).*s2-(k2*c1+k3*c1+k1*(c2+c3)).*s-k1*(k2+k3);
N=(m1.*s2+c1.*s+k1).*B + (c1.*s+k1).*C; 
Dstiff=N./B;

figure(1), hold on
plot (f, 20*log10(abs(Dstiff)),...
    'b-.', 'LineWidth', 3),
set(gca, 'XScale', 'log'),
xlim([0 3e4]), ylim([90 200]), grid on

E=(k2*h/S)


%%
%Plotto le masse dinamiche:
% m1=15.9; %massa piastra acciaio [kg];
% m2=21.768; %massa din. slab lab strade [kg];
% m3=1; %massa din. pavimento [kg];
omega=2*pi.*f;
M1=(m1+m2+m3).*omega.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M1)), 'b-.', 'LineWidth', 1),
set(gca, 'XScale', 'log'),

M2=(m1+m2).*omega.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M2)), 'g-.', 'LineWidth', 1),
set(gca, 'XScale', 'log'),

M3=m1.*omega.^2;
figure(1), hold on, 
plot (f, 20*log10(abs(M3)), 'k-.', 'LineWidth', 1),
set(gca, 'XScale', 'log'),
return

