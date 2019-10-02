%Martello strumentato 

set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables
load dati.mat

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
campione='Slab, resina 1gg';
accelerometro=1;
punta='m'; %m=metallica; p=plastica; g=gomma.
piastra='grande';
martellamento='auto';

 x = slab_piccola_m_res_1gg (:,1); % Force [N] 
 y = slab_piccola_m_res_1gg (:,accelerometro+2); % Accelerazione [m/s^2] 
%x = sito_piastra_piccola_lato_metallo_mart; % Force [N] 
%y = sito_piastra_piccola_lato_metallo_acc; % Accelerazione [m/s^2] 


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
soglia=10;               % Soglia dei picchi;
delay=round(0.1*fs);    % Sample da saltare una volta superata la soglia
inizio=10*fs;          % Punto di inizio della ricerca dei picchi;
fine=round(1*size(x));           % Punto di fine della ricerca dei picchi
% Parametri di filtro
bandwidth=400;            % Larghezza di banda richiesta al singolo colpo
% Dimensioni dei campioni
L_pre=round(1/2000*fs); % Lunghezza della parte prima del picco
L_coda=round(1*fs);   % Lunghezza della coda dei segnali
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
ascissamax=2000;        % Frequenza massima da plottare nei grafici
misura=['Campione ',num2str(campione),', ',martellamento,', punta ',punta,', piastra ',piastra,',Hann, PSDvsFFT, ',num2str(bandwidth),' Hz'];

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calibrazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
x=  div_F*x;     % Force [N] 
y=g*div_A*(-y);

L = length(x);
dt=1/fs; time=1000.*[0:dt:L/fs-dt];
figure(1)
subplot(2,1,1), hold on, plot (time, x)
subplot(2,1,2), hold on, plot (time, y), 
hold off

%<<<<<<<<<<<<<<<<<<<<
% Ricerca dei PICCHI
%<<<<<<<<<<<<<<<<<<<<

[picchi,n_picchi] = trovacolpi(x, soglia, delay, inizio, fine);
n_picchi

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[F, pos] = creamatriceforza_noavg (x, picchi,n_picchi, L_pre, L_coda, filt_doppi, fs);
[A] = creamatriceaccelerazione (y, pos, L_pre, L_coda, fs);
picchi_sel1=length(pos)

% Stampo dei picchi selezionati
x_sel = reshape(F,[],1);
y_sel = reshape(A,[],1);
l = length(x_sel);
dt=1/fs; time=1000*(0:dt:l/fs-dt);
figure(2) 
subplot(2,1,1), hold on, plot(time, x_sel)
subplot(2,1,2), hold on, plot(time, y_sel), 
hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<
%Faccio un calcolo di F_filt per ottenere L_win 
[F_filt, L_win] = finestra_forza (F, window_F, fs);

%Finestro sia Accelerazione che Fora utilizzando finestra_accelerazione 5
%sulla base di L_win
[A_filt] = finestra_accelerazione (A, window_A, L_win, fs);
[F_filt] = finestra_accelerazione (F, window_A, L_win, fs);
    
%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle FRFs
%<<<<<<<<<<<<<<<<<<<<
L = length(F_filt(:,1));
[PSD_F, f]= periodogram(F_filt, [], L, fs); %PSD Forza [N^2]
[PSD_A, f]= periodogram(A_filt, [], L, fs); %PSD Accelerazione [g^2]

%<<<<<<<<<<<<<<<<<<
% Filtraggio Banda
%<<<<<<<<<<<<<<<<<<
[r,c]=size(PSD_F);
tagli=[];
scarti=0;
for jj=1:(c-scarti)
    f0=find(f>ascissamin,1);
    fmax=find(PSD_F(f0:end, jj-scarti)<((PSD_F(f0, jj-scarti)/10)),1); 
    fmax=f(fmax+f0);
         if  fmax<bandwidth
            PSD_F(:,jj-scarti)=[];
            tagli=[tagli; jj];
            scarti=scarti+1;
         end
end
picchi_sel2= picchi_sel1 - scarti
PSD_A(:,tagli)=[];

% Matrici F e An filtrate con il bandwidth
F_filt2 = F_filt ;
A_filt2 = A_filt;
F_filt2(:,tagli)=[];
A_filt2(:,tagli)=[];

%Calcolo e memorizzo Spettri in fft di F, A, V e D
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% calcolo la media degli spettri
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%PSD
PSD_Fav = mean(sqrt(PSD_F), 2);
PSD_Aav = mean(sqrt(PSD_A), 2);
PSD_V1av = PSD_Aav./(2*pi*f); %velocità
PSD_D1av = PSD_V1av./(2*pi*f); % displacement

%FFT
FFT_Fav = mean( fft (F_filt2,  [], 1), 2); %FFT forza
FFT_Aav = mean( fft (A_filt2, [], 1), 2); %FFT accelerazione
FFT_V1av = FFT_Aav(1:L/2+1)./(1i*2*pi*f); %velocità
FFT_D1av=FFT_V1av(1:L/2+1)./(1i*2*pi*f); % displacement

%Ciclo for per plottare segnali e spettri (PSD):
dt=1/fs; time1=1000*(0:dt:L/fs-dt);
for j=1:(c-scarti)
        figure(101), grid on,
        sgtitle(misura)
        subplot(2,2,1), hold on, plot(time1, F_filt(:, j)), xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Force'), grid on, xlim([0 10])
        subplot(2,2,3), hold on, semilogx (f, 10*log10(PSD_F(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 N/Hz)'), title('PSD Force'), 
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        subplot(2,2,2), hold on, plot(time1, A_filt(:, j)), xlabel('Time [ms]'), ylabel('Amplitude [g]'), title('Acceleration'), grid on, xlim([0 10])
        subplot(2,2,4), hold on, semilogx (f, 10*log10(PSD_A(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration'), 
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
end
hold off

%Plot degli spettri medi in PSD
figure (101), subplot(2,2,3), hold on, plot (f, 20*log10(PSD_Fav), 'b-.', 'LineWidth', 3)
figure (101), subplot(2,2,4), hold on, plot (f, 20*log10(PSD_Aav), 'b', 'LineWidth', 3)
saveas (gcf, ['Segnali e spettri-C_',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.fig'])

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo della frequenza massima
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
PSD_Fav_dB = 20.*log10(PSD_Fav);
fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
fmax=f(fmax(1)+f0);
%plot sulla PSD della forza
figure (101), subplot (2,2,3),hold on
xl=xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom'; hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% COERENZA usando Forza / Accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
F_filtall = reshape(F_filt2, [],1);
A_filtall = reshape(A_filt2, [],1);
[r,c]=size(F_filt2);
[Cxy1,f] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],L,fs);
save (['Coherence, misura-C',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Cxy1');

clear F_filt2
clear A_filt2

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo DYNAMIC MASS Mechanical Impedance Dynamic Stiffness
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Dynamic Mass
DMASS1_av = PSD_Fav./PSD_Aav; %Modulo della Dynamic Mass;
% save (['DMASS1_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwidth),'Hz','.mat'], 'DMASS1_av');
DMASS1_av_ph = angle(FFT_Fav./FFT_Aav); %trovo la fase media usando le FFT;
% %<<<<<<<<<<<<<<<<<<<
% % Plot Dynamic Mass
% %<<<<<<<<<<<<<<<<<<<
% figure(104), 
% sgtitle(misura)
% subplot(3,1,1), plot(f,Cxy1), set(gca, 'XScale', 'log'), 
% title('Magnitude-Squared Coherence')
% xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
% grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 
% figure (104), subplot(3,1,2), hold on, semilogx (f, 20*log10(DMASS1_av), 'LineWidth', 3),
% % k=12.5e6;
% % semilogx (f, 20*log10(k./((2*pi*f').*(2*pi*f'))), 'LineWidth', 3),
% set(gca, 'XScale', 'log'), 
% xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Mass| (dB ref 1 kg)'), title(['Dynamic Mass (Force/Acceleration) Amplitude']), 
% xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']), grid on, xlim([ascissamin ascissamax])
% figure (104), subplot(3,1,3), 
% hold on,
% plot (f, 180.*DMASS1_av_ph(1:L/2+1)./pi, 'LineWidth', 3),
% set(gca, 'XScale', 'log')
% xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Mass (Force/Acceleration) Phase']),
% xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']);
% grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
% % saveas (gcf, ['Coerenza e dmass-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwidth),'Hz','.fig'])

% Mechanical Impedance
MI1_av = PSD_Fav./PSD_V1av; %Modulo dell'Impadenza meccanica
% save (['MI_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwidth),'Hz','.mat'], 'MI1_av');
MI1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_V1av); %trovo la fase media usando le FFT;
%<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Mechanical Impedence
%<<<<<<<<<<<<<<<<<<<<<<<<<<<
% figure(105), 
% sgtitle(misura)
% subplot(3,1,1), plot(f,Cxy1), set(gca, 'XScale', 'log'), 
% title('Magnitude-Squared Coherence')
% xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
% grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 
% figure (105), subplot(3,1,2), hold on, plot (f, 20.*log10(MI1_av), 'LineWidth', 3), 
% %semilogx (f, 20*log10(MI_av), 'LineWidth', 3),
% set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylim([0 130])
% xlabel('log(Frequency) [Hz]'), ylabel('20 log |Mech. Impedance| (dB ref 1 N s/m]'), title(['Mechanical Impedance (Force/Velocity) Amplitude']), 
% xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']); grid on, xlim([ascissamin ascissamax])
% figure (105), subplot(3,1,3), 
% hold on, plot (f, 180.*MI1_av_ph./pi, 'LineWidth', 3),
% set(gca, 'XScale', 'log')
% xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Mechanical Impedance (Force/Velocity) Phase']),
% xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']);
% grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
% % saveas (gcf, ['Coerenza e MI-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwidth),'Hz','.fig'])


%_____________________________________________________
% Dynamic Stiffness
Dstiff1_av = PSD_Fav./PSD_D1av; %Modulo della Dynamic Stiffness
save (['Dstiffness_av, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av');
Dstiff1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_D1av); %trovo la fase media usando le FFT;
save (['Dstiffness_av_ph, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av_ph');

%<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dynamic Stiffness
%<<<<<<<<<<<<<<<<<<<<<<<<
figure(106), 
sgtitle(misura)
subplot(3,1,1), plot(f,Cxy1), set(gca, 'XScale', 'log'), 
title('Magnitude-Squared Coherence')
xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 
figure (106), subplot(3,1,2), hold on, semilogx (f, 20*log10(Dstiff1_av), 'LineWidth', 3),
%semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylim([100 200])
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'), title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']); grid on, xlim([ascissamin ascissamax])
figure (106), subplot(3,1,3), 
hold on,
plot (f, 180.*Dstiff1_av_ph./pi, 'LineWidth', 3),
set(gca, 'XScale', 'log')
xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Stiffness (Force/Displacement) Phase']),
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']);
grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
saveas (gcf, ['Coerenza e Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,', ',num2str(bandwidth),'Hz','.fig'])                    

%<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot singolo D-Stiffness
%<<<<<<<<<<<<<<<<<<<<<<<<<<

figure (107),
semilogx (f, 20*log10(Dstiff1_av), 'LineWidth', 3),
%semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'), title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']);
grid on, xlim([20 ascissamax]),ylim([120 190])


%%
%Derivo il m,odulo elsstico dinamico della pavimentazione ricavando il massimo della Dstiff 
% e usando come spessore della pavimentazione 4 cm e l'area di carico circolare con
% raggio pari a 7.5 mm (piastina);
c0=1;
c1=2*pi*0.1;
c2=2*pi*0.2;
c3=2*pi*0.3;
c4=2*pi*0.6;
c5=2*pi*0.7;
c=[c0 c1 c2 c3 c4 c5];

C0m=1;
C1m=1;
C2m=0.8872;
C3m=1.2596;
C4m=1.0883;
C5m=1.0508;
C=[C0m C1m C2m C3m C4m C5m];

% %kmax=max(Dstiff1_av (300:600,1))/(C(accelerometro+1)*c(accelerometro+1));
% kmax=max(Dstiff1_av (300:1000,1));
% E= kmax*(0.04) / (pi*0.0075^2) /1e6



%%
%Risultato 10.09.19:
% con piastra di carico piccola
% 8 picchi selezionati e bandwidth = 434 Hz
% E = ((10^((164.2)/20))*(0.04)/(0.05*0.056))/1e6 = 2317 MPa
% Dmass=10^(16.76/20)=6.9 kg





%__________________________________________________________
%Risultati 09.09.19:
%spessore pavim.=4 cm; raggio piastrina=7.5mm;
%Punto 1
%K=148.9 dB(N/m)
%f0=416.8 Hz
%E1.1=6306 MPa;
%_________________
%K= 149.8 dB(N/m)
%f0= 424.8 Hz
%E1.2= 6995 MPa;
%_________________
%K=  148.5 dB(N/m)
%f0=  411.8 Hz
%E1.3= 6023 MPa;
%__________________________media= 6441 MPa

%Punto 2
%K= 155.1 dB(N/m)
%f0= 472.8 Hz
%E2.1 =  12876 MPa;
%_________________
%K=  145.6 dB(N/m)
%f0=  418.8 Hz
%E2.2 =  4313  MPa;
%_________________
%K=  149.1 dB(N/m)
%f0=  437.8 Hz
%E2.3 =  6453 MPa;
%__________________________media=7881 MPa

%Punto 3
%K= 146  dB(N/m)
%f0=  452.8 Hz
%E3.1 =  4516 MPa;
%_________________
%K= 147.7  dB(N/m)
%f0=  408.8 Hz
%E3.2 =  5493 MPa;
%_________________
%K= 149.5 dB(N/m)
%f0= 398.8 Hz
%E3.3 =  6757 MPa;
%_______________________media= 5589 MPa







%%
%vecchi risultati:
%1.1
E0 = 6.0159e+03;
E1=1.9101e+05;
E2=1.6024e+05;
E3=1.4775e+05;
E4=  2.5726e+05;
E5= 2.5017e+05;

%1.2
E0 =  6.6426e+03;

%1.3
E0 = 3.7155e+03;

%2.1
 E0=9.4884e+03;
%2.2
E0 =5.6057e+03
%2.3
E = 5.0330e+03


%%
E= ((10^((150.6)/20))*(0.04)/(pi*0.0075^2))/1e6

%prove del 21.08.19, modulo elastico della piastra di carico 2 (su spugna),
%la piastra di spessore 9 mm, considerando un K=177 N/m, 
%e un'area di carico di raggio 3 mm, è risultato 225 GPa (valore
%di letteratura = 200 GPa per steel).

%prove del 21.08.19, modulo elastico della piastra di carico 2 (su pavimento),
%la piastra di spessore 9 mm, K=23 N/m

%prove del 21.08.19, modulo elastico della piastra di carico 2 (su pavimento),
% K=187 N/m

%prove del 21.08.19, modulo elastico della piastra di carico 2 (su pavimento),
% K=173 N/m


%%
% clc
% close all
% clear all

Dstiff_av_mod_a1=abs(Dstiff_av1);
Dstiff_av_mod_a2=abs(Dstiff_av2);
figure (108), subplot(2,1,1), hold on, loglog (f, Dstiff_av_mod_a2, 'LineWidth', 3), set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
figure (108), subplot(2,1,1), hold on, loglog (f, Dstiff_av_mod_a1, 'LineWidth', 3), set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
legend('Accelerometro 2 (picchiato male)','Accelerometro 1')
hold off

figure (107),hold on, loglog(f,Dstiff_av_mod_a1./Dstiff_av_mod_a2),
% set(gca, 'XScale', 'log'), 
% set(gca, 'YScale', 'log')
xlim([ascissamin ascissamax]),hold off,

xlabel('log(Frequency) [Hz]'), ylabel('Amplitude [N / m]'), title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']); grid on, xlim([ascissamin ascissamax])
legend('Media sulle PSD')
figure (106), subplot(2,1,2), 
hold on,
plot (f, 180.*Dstiff_av_ph(1:L/2+1)./pi, 'LineWidth', 3),
set(gca, 'XScale', 'log')
xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Stiffness (Force/Displacement) Phase']),
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']);
grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
%saveas (gcf, ['Coerenza e Dstiff-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.fig'])                    


%%
close all
%DSTIFF=[Dstiff_1_1, Dstiff_1_2, Dstiff_1_3, Dstiff_2_1, Dstiff_2_2,Dstiff_2_3, Dstiff_3_1, Dstiff_3_2, Dstiff_3_3]; 
D=mean(DSTIFF, 2);
for i=1:9
figure(1), hold on, plot(f, 20*log10(DSTIFF(:,i)))
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
grid on, xlim([100 2000])
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'), title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 

end
plot(20*log10(D), 'r-.', 'LineWidth', 3)
lgd=legend ({'1.1', '1.2', '1.3', '2.1', '2.2', '2.3', '3.1', '3.2', '3.3', 'media', 'piastra piccola'}, 'Location', 'west')
lgd.FontSize = 16;
hold off

%%
% Confronti tra ripetizioni
%clear variables

figure
hold on 
titolo= 'Confronto biadesivo resina nel tempo, lato';
sgtitle(titolo)
% colori: #0072BD #D95319 #EDB120 #7E2F8E #77AC30

% %confront pavimentazioni
% %Via cocchi
% semilogx (f, 20*log10(Dstiff1_av_G_1_diretta),'-','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_2_diretta),'-','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_3_diretta),'-','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_4_diretta),'-','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_5_diretta),'-','Color','#0072BD','LineWidth', 2),
% 
% semilogx (f, 20*log10(Dstiff1_av_G_1_trasmessa),'--','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_2_trasmessa),'--','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_3_trasmessa),'--','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_4_trasmessa),'--','Color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_5_trasmessa),'--','Color','#0072BD','LineWidth', 2),
% %via battelli3 1
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_1_diretta),'-','Color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_2_diretta),'-','Color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_3_diretta),'-','Color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_4_diretta),'-','Color','#D95319','LineWidth', 2),
% 
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_1_trasmessa),'--','Color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_2_trasmessa),'--','Color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_3_trasmessa),'--','Color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli1_4_trasmessa),'--','Color','#D95319','LineWidth', 2),
% 
% %via battelli3 2
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_1_diretta),'-','Color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_2_diretta),'-','Color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_3_diretta),'-','Color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_4_diretta),'-','Color','#EDB120','LineWidth', 2),
% 
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_1_trasmessa),'--','Color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_2_trasmessa),'--','Color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_3_trasmessa),'--','Color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelli2_4_trasmessa),'--','Color','#EDB120','LineWidth', 2),
% 
% %via battelli3 parck
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_1_diretta),'-','Color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_2_diretta),'-','Color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_3_diretta),'-','Color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_4_diretta),'-','Color','#7E2F8E','LineWidth', 2),
% 
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_1_trasmessa),'--','Color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_2_trasmessa),'--','Color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_3_trasmessa),'--','Color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_battelliparck_4_trasmessa),'--','Color','#7E2F8E','LineWidth', 2),
% legend('Via Cocchi "diretta"','','','','','Via Cocchi "trasmessa"','','','','','Via battelli 3,1 "diretta"','','','','Via battelli 3,1 "trasmessa"','','','','Via battelli 3,2 "diretta"','','','','Via battelli 3,2 "trasmessa"','','','','Via battelli parcking "diretta"','','','','Via battelli parcking "trasmessa"')


% confronto Resina
semilogx (f, 20*log10(Dstiff1_av_P_bi),'--','LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_5 ), 'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_10), 'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_20), 'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_40), 'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_80), 'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_160),'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_300),'LineWidth', 2),
semilogx (f, 20*log10(Dstiff1_av_P_res_1gg),'LineWidth', 2),

legend('Piccola biad','Piccola Resina 5 min','Piccola Resina 10 min','Piccola Resina 20 min','Piccola Resina 40 min','Piccola Resina 80 min','Piccola Resina 160 min','Piccola Resina 300min','Piccola Resina 1 gg')


% % confronto Pattex Millechiodi
% semilogx (f, 20*log10(Dstiff1_av_P_bi),'LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_P_at),'LineWidth',2),
% semilogx (f, 20*log10(Dstiff1_av_P_patt),'LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_bi),'--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_at),'--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_patt),'--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_G_patt2),':','LineWidth', 2),
 
% semilogx (f, 20*log10(Dstiff1_av_1G),'b','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_2G),'b','LineWidth',2),
% semilogx (f, 20*log10(Dstiff1_av_3G),'b','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_1P),'b--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_2P),'b--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_3P),'b--','LineWidth', 2),

% semilogx (f, 20*log10(Dstiff1_av_1Gn),'k','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_2Gn),'k','LineWidth',2),
% semilogx (f, 20*log10(Dstiff1_av_3Gn),'k','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_1Pn),'k--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_2Pn),'k--','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff1_av_3Pn),'k--','LineWidth', 2),

%semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'),
%ylim([100 200])
xlabel('log(Frequency) [Hz]')
ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
%title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
grid on,
xlim([20 2000]), %ylim([130 170])

%legend('Piccola biad','Piccola Attak','Piccola pattex','Grande biadesivo','Grande attak','Grande pattex','Grande pattex 2gg')

%legend('Grande Rip. 1','Grande Rip. 2','Grande Rip. 3','Piccola Rip. 1','Piccola Rip. 2','Piccola Rip. 3','Grande Rip. 1 no av','Grande Rip. 2 no av','Grande Rip. 3 no av','Piccola Rip. 1 no av','Piccola Rip. 2 no av','Piccola Rip. 3 no av')
saveas (gcf, [titolo,'.fig'])




