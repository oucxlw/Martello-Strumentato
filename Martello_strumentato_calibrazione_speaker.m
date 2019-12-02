%Martello strumentato 

set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables
load segnali_accelerazione.mat
scarti=0;
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
campione='Calibrazione speaker';
accelerometro=0;
punta=''; %m=metallica; p=plastica; g=gomma.
piastra='speaker';
martellamento='auto';

 %x = cal_biad_1 (:,1); % Force [N] 
 y0 = cal_biad_0 (:,accelerometro+2); % Accelerazione [m/s^2]
 y1 = cal_biad_1 (:,accelerometro+3);
 y4 = cal_biad_4 (:,accelerometro+6);
 
%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri di controllo
%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri fisici
% sensit_F=1/1000;      % Sensibilità Martello strumentato [V/N]
% sensit_A1=1/1000;%95.8/1000;    % Sensibilità Accelerometro1 [V/g]
% sensit_A2=1/1000;%0.03/1000;   % Sensibilità Accelerometro2 [V/g]
g=9.81;                 % Accelerazione di gravità [m/s^2]
div_F=2000;             % Divisione applicata alla Forza prima della scrittura sul file WAVE
div_A=500;              % Divisione applicata alla Accelerazione prima della scrittura sul file WAVE

% Parametri di ricerca
fs=52100;               % Freq. di campionamento (Hz);
soglia0=0.04;
soglia1=0.004;
soglia4=0.004; % Soglia dei picchi;
delay=round(1*fs);    % Sample da saltare una volta superata la soglia

inizio0=470000;          % Punto di inizio della ricerca dei picchi;
fine0=round(1*size(y0));           % Punto di fine della ricerca dei picchi
inizio1=488000;
fine1=round(1*size(y1));
inizio4=693000;
fine4=round(1*size(y4));

% Parametri di filtro
bandwidth=400;            % Larghezza di banda richiesta al singolo colpo
% Dimensioni dei campioni
L_pre=round(1/100*fs); % Lunghezza della parte prima del picco
L_coda=round(30*fs);   % Lunghezza della coda dei segnali
% Filtraggio doppi colpi
filt_doppi=0;           % Se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi
% Normalizzazione colpi
norm=0;                 % Se norm=1 i colpi vengono normalizzati
% Finestratura
window_F=5;         % Tipo di finestratura da applicare
window_A=5;
% 0 = nessuna finestratura
% 1 = finestratura quadrata
% 2 = hamming (da verificare)
% 3 = blackmanharris ritardata
% 4 = blackmanharris anticipata
% 5 = Hann Window
% Plotting
ascissamin=50;         % Frequenza minima da plottare nei grafici 
ascissamax=1000;        % Frequenza massima da plottare nei grafici
misura=['Campione ',num2str(campione),', ',martellamento,', punta ',punta,', piastra ',piastra,',Hann, PSDvsFFT, ',num2str(bandwidth),' Hz'];

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calibrazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%x=  div_F*x;     % Force [N] 
y0=g*div_A*(-y0);
y1=g*div_A*(-y1);
y4=g*div_A*(-y4);

L0 = length(y0);
dt=1/fs; time=1000.*[0:dt:L0/fs-dt];
figure(1)
subplot(3,1,1), hold on, plot (y0),
subplot(3,1,2), hold on, plot (y1),
subplot(3,1,3), hold on, plot (y4),
hold off




%<<<<<<<<<<<<<<<<<<<<
% Ricerca dei PICCHI
%<<<<<<<<<<<<<<<<<<<<

[picchi0,n_picchi0] = trovacolpi(y0, soglia0, delay, inizio0, fine0);
[picchi1,n_picchi1] = trovacolpi(y1, soglia1, delay, inizio1, fine1);
[picchi4,n_picchi4] = trovacolpi(y4, soglia4, delay, inizio4, fine4);


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[A0, pos] = creamatriceforza_noavg (y0, picchi0, 1, L_pre, L_coda, filt_doppi, fs);
[A1, pos] = creamatriceforza_noavg (y1, picchi1, 1, L_pre, L_coda, filt_doppi, fs);
[A4, pos] = creamatriceforza_noavg (y4, picchi4, 1, L_pre, L_coda, filt_doppi, fs);

%[A] = creamatriceaccelerazione (y, pos, L_pre, L_coda, fs);
picchi_sel1=length(pos)

% Stampo dei picchi selezionati
%x_sel = reshape(F,[],1);
y0_sel = reshape(A0,[],1);
y1_sel = reshape(A1,[],1);
y4_sel = reshape(A4,[],1);

l = length(y0_sel);
dt=1/fs; time=1000*(0:dt:l/fs-dt);
figure(2)
%subplot(2,1,1), hold on, plot(time, x_sel)
subplot(3,1,1), hold on, plot(time, y0_sel),
subplot(3,1,2), hold on, plot(time, y1_sel), 
subplot(3,1,3), hold on, plot(time, y4_sel), 

hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<
%Faccio un calcolo di F_filt per ottenere L_win 
%[F_filt, L_win] = finestra_forza (F, window_F, fs);

%Finestro sia Accelerazione che Fora utilizzando finestra_accelerazione 5
%sulla base di L_win
[A0_filt] = A0;%finestra_accelerazione (A, window_A, L_win, fs);
[A1_filt] = A1;
[A4_filt] = A4;
%[F_filt] = finestra_accelerazione (F, window_A, L_win, fs);
    
%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle FRFs
%<<<<<<<<<<<<<<<<<<<<
L0 = length(A0_filt(:,1));
%[PSD_F, f]= periodogram(F_filt, [], L, fs); %PSD Forza [N^2]
L1 = length(A1_filt(:,1));
L4 = length(A4_filt(:,1));

[PSD_A0, f]= periodogram(A0_filt, [], L0, fs); %PSD Accelerazione [g^2]
[PSD_A1, f]= periodogram(A1_filt, [], L1, fs); %PSD Accelerazione [g^2]
[PSD_A4, f]= periodogram(A4_filt, [], L4, fs); %PSD Accelerazione [g^2]


A_filt2 = A0_filt;
% F_filt2(:,tagli)=[];
% A_filt2(:,tagli)=[];

%Calcolo e memorizzo Spettri in fft di F, A, V e D
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% calcolo la media degli spettri
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%PSD

PSD_A0av = mean(sqrt(PSD_A0), 2);
PSD_A1av = mean(sqrt(PSD_A1), 2);
PSD_A4av = mean(sqrt(PSD_A4), 2);

%FFT
%FFT_Fav = mean( fft (F_filt2,  [], 1), 2); %FFT forza
FFT_Aav = mean( fft (A_filt2, [], 1), 2); %FFT accelerazione
FFT_V1av = FFT_Aav(1:L0/2+1)./(1i*2*pi*f); %velocità
FFT_D1av=FFT_V1av(1:L0/2+1)./(1i*2*pi*f); % displacement

%Ciclo for per plottare segnali e spettri (PSD):
dt=1/fs; time1=1000*(0:dt:L0/fs-dt);

% for j=1:(c-scarti)
%         figure(101), grid on,
%         sgtitle(misura)
%         %subplot(2,2,1), hold on, plot(time1, F_filt(:, j)), xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Force'), grid on, xlim([0 10])
%         %subplot(2,2,3), hold on, semilogx (f, 10*log10(PSD_F(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 N/Hz)'), title('PSD Force'), 
%         grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
%         subplot(1,2,1), hold on, plot(time1, A0_filt(:, j)), xlabel('Time [ms]'), ylabel('Amplitude [g]'), title('Acceleration'), grid on, xlim([0 10])
%         subplot(1,2,2), hold on, semilogx (f, 10*log10(PSD_A0(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration'), 
%         grid on, set(gca, 'XScale', 'log'), xlim([50 1000])
% end
% hold off

figure(101), grid on,
sgtitle(misura)
%subplot(2,2,1), hold on, plot(time1, F_filt(:, j)), xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Force'), grid on, xlim([0 10])
%subplot(2,2,3), hold on, semilogx (f, 10*log10(PSD_F(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 N/Hz)'), title('PSD Force'), 
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

subplot(3,1,1), hold on, grid on
semilogx (f, 10*log10(PSD_A0av)),
xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration 0'), 
set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax]),ylim([-50 -10]),

subplot(3,1,2), hold on, grid on
semilogx (f, 10*log10(PSD_A1av)),
xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration 1'),
set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax]),ylim([-50 -10]),

subplot(3,1,3), hold on, grid on
semilogx (f, 10*log10(PSD_A4av)),
xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration 4'), 
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax]),ylim([-50 -10]),

C01 = sqrt(PSD_A0av./PSD_A1av);
C04 = sqrt(PSD_A0av./PSD_A4av);
C14 = sqrt(PSD_A1av./PSD_A4av);


figure
hold on
plot (f,C01)
plot (f,C04)
plot (f,C14)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax]), ylim([0 2])
legend('C01','C04','C14')

%Plot degli spettri medi in PSD
%figure (101), subplot(2,2,3), hold on, plot (f, 20*log10(PSD_Fav), 'b-.', 'LineWidth', 3)
%figure (101), subplot(2,2,4), hold on, plot (f, 20*log10(PSD_Aav), 'b', 'LineWidth', 3)
saveas (gcf, ['Segnali e spettri-C_',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.fig'])
return
%%

