%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Modello matematico per i campioni di piccole dimensioni
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

clear all
%close all
load f %vettore f

f=[f;f(1:end)+f(end)+f(2)] % Allungo f fino a fs
w=2*pi*f; % Creo vettore omega
n=3; %ordine del sistema
m=ones(n,1);
k=10^9*ones(n,1);
z=ones(n,1);

% data=[1,  50*10^4,   0.03,  20e+3,  3,  16*pi];
% m=[data(1),data(3),data(5)];
% k=[data(2),data(4),data(6)];

% c=z/100;
for ii=0.1%:0.1:0.9

m1=1.4;
ripartizione=0.3;
m=[m1*(1-ripartizione), m1*ripartizione,0.038+ 7];
k=[0.35*2.8*10^7, 0.15*2.8*10^6, 10^1]; %k=[5*10^10 2.8*10^8 10^2];
c=k.*[0.00000000001, 0.0006, 1];
[K,I,M,B]= calculum(m,k,c,w);


fig1=figure(107);
hold on
plot(f,20*log10(abs(K)))
grid on,ylim([40 220]),
fig1.Children.XScale='log';

%clearvars K I M B
end

load Forza %In 'Forza' c'è F che è la forza nel tempo espressa in N
t=1:length(F);
t=t./52100;
i=30; % indice del segnale da plottare
figure
hold on
plot(t,F(:,i));

F_f=fft(F);
% plot(abs(F_f));

A=ifft(-F_f./M(1:end-1),length(F),'symmetric');
plot( t,A(:,i) )

save('Dati_simulazione.mat','F','A','K','I','M')
%save('pippo.mat','data','F',...)

load Forza_filtrata %In forza c'è F che è la forza nel tempo espressa in N
t=1:length(F_filt);
t=t./52100;
i=30;
figure
hold on
plot(t,F_filt(:,i));

F_f=fft(F_filt);
% plot(abs(F_f));

A=ifft(-F_f./M(1:end-1),length(F),'symmetric');
plot( t,A(:,i) )

% save('Dati_simulazione.mat','F','A','K','I','M')

%%

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

set(gca, 'XScale', 'log'),
set(gca, 'YScale', 'lin'),
