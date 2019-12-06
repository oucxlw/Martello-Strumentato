% Martello strumentato filtro intensita
%
set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables

%%
clc
close all
clear variables

load dati.mat

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
campione='Campione 1, piastra pesante 1, punta metallica, nessun collante';
accelerometro=0;
punta='M'; %m=metallica; p=plastica; g=gomma.
piastra='pesante 1';
martellamento='auto';

x = c1_a0_pesante1_met_terra_nonattaccata(:,1); % Force [N]
y = c1_a0_pesante1_met_terra_nonattaccata (:,2); % Accelerazione [m/s^2]


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
soglia=10;              % Soglia dei picchi;
delay=round(0.5*fs);    % Sample da saltare una volta superata la soglia
inizio=1*fs;            % Punto di inizio della ricerca dei picchi;
fine=round(0.95*size(x));           % Punto di fine della ricerca dei picchi
% Parametri di filtro
bandwidth=0;            % Larghezza di banda richiesta al singolo colpo
% Dimensioni dei campioni
<<<<<<< Updated upstream
L_pre=round(1/2000*fs); % Lunghezza della parte prima del picco
L_coda=round(1*fs);     % Lunghezza della coda dei segnali
=======
L_pre=round((2^15/2)); % Lunghezza della parte prima del picco
L_coda=round(2^15-L_pre-1);     % Lunghezza della coda dei segnali
>>>>>>> Stashed changes
% Filtraggio doppi colpi
filt_doppi=0;           % Se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi
% Normalizzazione colpi
norm=0;                 % Se norm=1 i colpi vengono normalizzati
% Finestratura
window_F=5;             % Tipo di finestratura da applicare
window_A=5;
% 0 = nessuna finestratura
% 1 = finestratura quadrata
% 2 = hamming (da verificare)
% 3 = blackmanharris ritardata
% 4 = blackmanharris anticipata
% 5 = Hann Window
% Plotting
<<<<<<< Updated upstream
ascissamin=1;         % Frequenza minima da plottare nei grafici 
ascissamax=20000;       % Frequenza massima da plottare nei grafici
=======
ascissamin=1;         % Frequenza minima da plottare nei grafici
ascissamax=5000;       % Frequenza massima da plottare nei grafici
>>>>>>> Stashed changes
misura=['Campione ',num2str(campione),', ',martellamento,', punta ',punta,', piastra ',piastra,',Hann, PSDvsFFT, ',num2str(bandwidth),' Hz'];

%colore=['#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'; '#4DBEEE'; '#A2142F';'#6495ed'; '#b8860b'; '#a9a9a9'; '#cd5c5c'; '#20b2aa'; '#4DBEEE'; '#A2142F'];
colore=[
    '.000, .447, .741';
    '.850, .325, .098';
    '.929, .694, .125';
    '.494, .184, .556';
    '.466, .674, .188';
    '.301, .745, .933';
    '.635, .078, .184';
    '.762, .762, .762';
    '.999, .682, .788';
    '.000, .502, .502';
    '.000, .000, .627';
    '.430, .135, .078';
    '.559, .354, .075';
    '.424, .124, .526';
    '.446, .644, .148';
    '.301, .705, .903';
    '.615, .018, .114'];
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calibrazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
x=  div_F*x;     % Force [N]
y=g*div_A*(-y);

L = length(x);
dt=1/fs; time=[0:dt:L/fs-dt];
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

[r,c]=size(F);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% F=[]; A=[];
% load Dati_simulazione
% A=real(A);

<<<<<<< Updated upstream
%Faccio un calcolo di F_filt per ottenere L_win 
[~, L_win] = finestra_forza (F, window_F, fs);

%Finestro sia Accelerazione che Forza utilizzando finestra_accelerazione 5
%sulla base di L_win
[A_filt] = finestra_accelerazione (A, window_A, L_win, fs);
[F_filt] = finestra_accelerazione (F, window_A, L_win, fs);

%<<<<<<<<<<<<<<<<<<<<<<<<<
% Abbattimento della coda
%<<<<<<<<<<<<<<<<<<<<<<<<<
divid=20;
for ii=1:CC
    L=L_win(ii); %PER CIASCUNA COLONNA DI f PRENDO LA LUNGHEZZA APPROPRIATA NEL VETTORE LUNGHEZZE L_win
    F_filt(:,ii)=[F_filt(1:(L-round(L/divid)-1),ii); movavg(F_filt((L-round(L/divid)):end,ii),'linear',round(L/7))];
end
=======
% %Faccio un calcolo di F_filt per ottenere L_win
% [~, L_win] = finestra_forza (F, window_F, fs);
% 
% %Non applico una vera finestratura.
% [A_filt] = finestra_accelerazione (A, window_A, L_win, fs);
% [F_filt] = finestra_accelerazione (F, window_A, L_win, fs);

% win=hann(length(A(:,1)));%[0;ones(r-2,1);0]; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

for i=15
M_win=round(2^i);
win=[zeros((r-M_win)/2,1);hann(M_win);zeros((r-M_win)/2,1)];
E_win=sum(win.^2)/length(win);

%finestratura
A_filt=A.*win;
F_filt=F.*win;
%trasformata (non normalizzata)
FFT_F=fft(F_filt);
FFT_A=fft(A_filt);
f_fft=0:(r-1);
f_fft=f_fft'/(r-1)*fs;


PSD_F_fft = 2*abs(FFT_F).^2./length (FFT_F(:,1))/fs/E_win; %PSD di F stimata tramite fft
PSD_A_fft = 2*abs(FFT_A).^2./length (FFT_A(:,1))/fs/E_win;

% %<<<<<<<<<<<<<<<<<<<<<<<<<
% % Abbattimento della coda
% %<<<<<<<<<<<<<<<<<<<<<<<<<
% divid=20;
% for ii=1:c
%     L=L_win(ii); %PER CIASCUNA COLONNA DI f PRENDO LA LUNGHEZZA APPROPRIATA NEL VETTORE LUNGHEZZE L_win
%     F_filt(:,ii)=[F_filt(1:(L-round(L/divid)-1),ii); movavg(F_filt((L-round(L/divid)):end,ii),'linear',round(L/7))];
% end
>>>>>>> Stashed changes

%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle PSDs
%<<<<<<<<<<<<<<<<<<<<
% PSD_win=ones(r,1);
[PSD_F, f]= periodogram(F_filt, win, r, fs); %PSD Forza [N^2]
[PSD_A, f]= periodogram(A_filt, win, r, fs); %PSD Accelerazione [g^2]
save ('f.mat', 'f');

%<<<<<<<<<<<<<<<<<<<<<<
% Filtraggio per Banda
%<<<<<<<<<<<<<<<<<<<<<<
[R,C]=size(PSD_F);
tagli=[];
scarti=0;
for jj=1:(C-scarti)
    f0=find(f>ascissamin,1);
    fmax=find(PSD_F(f0:end, jj-scarti)<((PSD_F(f0, jj-scarti)/10)),1);
    fmax=f(fmax+f0);
    if  fmax<bandwidth
        PSD_F(:,jj-scarti)=[];
        tagli=[tagli; jj];
        scarti=scarti+1;
    end
end
picchi_sel2 = picchi_sel1 - scarti
PSD_A(:,tagli)=[];

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi statistica pre filtraggio
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
bin=round(sqrt(picchi_sel2))+1;

[Y,E] = discretize(sqrt(max(PSD_F)),bin);
values=1:bin;

figure (3), hold on
histfit(sqrt(max(PSD_F)),bin);

%<<<<<<<<<<<<<<<<<<<<<<<
% Filtraggio Intensita'
%<<<<<<<<<<<<<<<<<<<<<<<
% [r,c]=size(PSD_F);
% tagli=[];
% scarti=0;
%
% for jj=1:(c-scarti)
%     if sqrt(max(PSD_F(:,jj-scarti)))<0.04 || sqrt(max(PSD_F(:,jj-scarti)))>0.16
%         PSD_F(:,jj-scarti)=[];
%         tagli=[tagli; jj];
%         scarti=scarti+1;
%     end
% end
% picchi_sel2 = picchi_sel2 - scarti
% PSD_A(:,tagli)=[];

%<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo Dstiff totale
%<<<<<<<<<<<<<<<<<<<<<<<

<<<<<<< Updated upstream
PSD_Fav = mean(sqrt(PSD_F), 2);
PSD_Aav = mean(sqrt(PSD_A), 2);
PSD_V1av = PSD_Aav./(1i*2*pi*f); %velocità
PSD_D1av = PSD_V1av./(1i*2*pi*f); % displacement

Dstiff=PSD_Fav./PSD_D1av;

PSD_F_fft = 2*abs(fft(F)).^2./length (F(:,1))/fs;
figure(110), hold on
subplot (2,1,1),plot(f,20*log10(PSD_F_fft(1:length(f),1)))
hold on,subplot (2,1,1),plot(f,20*log10(PSD_F(1:length(f),1)),'LineWidth',2)
,grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

PSD_A_fft = 2*abs(fft(A)).^2./length (A(:,1))/fs;
subplot (2,1,2),plot(f,20*log10(PSD_A_fft(1:length(f),1)))
hold on, subplot (2,1,2),plot(f,20*log10(PSD_A(1:length(f),1)),'LineWidth',2)
,grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
return
%<<<<<<<<<<<<<<<<<<<<<
% Calcolo tramite FFT
%<<<<<<<<<<<<<<<<<<<<<
f2=[f;f(2:end)+f(end)]; % Allungo f fino a fs
FFT_Fav = mean( fft(F),2);
FFT_Aav = mean( fft(A),2);
FFT_Vav = FFT_Aav./(1i*2*pi*f2); 
FFT_Dav = FFT_Aav./(1i*2*pi*f2).^2; 

FFT_Kav = FFT_Fav./FFT_Dav;
figure(109),hold on, 
subplot(3,1,2)
plot (f2,20*log10(abs(FFT_Kav)),'-','LineWidth',1)


=======
% FFT_Fav_pgram = sqrt(mean(PSD_F, 2));
% FFT_Aav_pgram = sqrt(mean(PSD_A, 2));
% FFT_Vav_pgram = FFT_Aav_pgram./(1i*2*pi*f); %velocità
% FFT_Dav_pgram = FFT_Vav_pgram./(1i*2*pi*f); % displacement
% FFT_Kav_pgram = FFT_Fav_pgram./FFT_Dav_pgram; %dynamic stiffness calcolata tramite periodogram

PSD_Fav_pgram = mean(PSD_F, 2);
PSD_Aav_pgram = mean(PSD_A, 2);
PSD_Vav_pgram = PSD_Aav_pgram./(1i*2*pi*f).^2; %velocità
PSD_Dav_pgram = PSD_Vav_pgram./(1i*2*pi*f).^2; % displacement
PSD_Kav_pgram = PSD_Fav_pgram./PSD_Dav_pgram; %dynamic stiffness calcolata tramite periodogram

% figure,hold on
% plot(F_filt(:,1))
% %<<<<<<<<<<<<<<<<<<<<<<<<<
% % Abbattimento della coda
% %<<<<<<<<<<<<<<<<<<<<<<<<<
% divid=20;
% for ii=1:C
%     L=L_win(ii); %PER CIASCUNA COLONNA DI f PRENDO LA LUNGHEZZA APPROPRIATA NEL VETTORE LUNGHEZZE L_win
%     F_filt(:,ii)=[F_filt(1:(L-round(L/divid)-1),ii); movavg(F_filt((L-round(L/divid)):end,ii),'linear',round(L/7))];
% end
% hold on
% plot(F_filt(:,1))

figure(110), hold on
subplot (2,1,1),hold on,
plot(f_fft,20*log10(PSD_F_fft(:,1)),'color',string(colore(i-7,:)))
plot(f,20*log10(PSD_F(1:length(f),1)),'-.','color',string(colore(i-7,:)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
legend('PSD tramite Fft','PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

subplot (2,1,2),hold on
plot(f_fft,20*log10(PSD_A_fft(:,1)),'color',string(colore(i-7,:)))
plot(f,20*log10(PSD_A(1:length(f),1)),'-.','color',string(colore(i-7,:)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

legend('PSD tramite Fft','PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

%<<<<<<<<<<<<<<<<<<<<<
% Calcolo tramite FFT
%<<<<<<<<<<<<<<<<<<<<<

FFT_Fav_fft = mean(FFT_F,2);
FFT_Aav_fft = mean(FFT_A,2);
FFT_Vav_fft = FFT_Aav_fft./(1i*2*pi*f_fft);
FFT_Dav_fft = FFT_Aav_fft./(1i*2*pi*f_fft).^2;
FFT_Kav_fft = FFT_Fav_fft./FFT_Dav_fft;
>>>>>>> Stashed changes

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% COERENZA usando Forza / Accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<< Updated upstream
F_filtall = reshape(F_filt, [],1);
A_filtall = reshape(A_filt, [],1);
[r,c]=size(F_filt);
[Cxy1,f] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],L,fs);
save ([num2str(round(indice)),'Coherence, misura-C',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Cxy1');
figure (109),hold on,
subplot(3,1,1)
semilogx(f,Cxy1)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

FFT_Kav_gate=FFT_Kav;
for i=1:length(f)
    if Cxy1(i)<0.5
        FFT_Kav_gate(i)=1;
        
    end
end
figure(109),hold on, 
subplot(3,1,2)
% semilogx (f2,20*log10(abs(FFT_Kav)),'LineWidth',1)
semilogx (f2,20*log10(abs(FFT_Kav_gate)),'r.','LineWidth',2)
grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]), ylim([100 200])
hold off
subplot(3,1,3)
semilogx (f2,180/pi*angle(FFT_Kav))
grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]);ylim([-180 180])





=======

F_filtall = reshape(F_filt, [],1);
A_filtall = reshape(A_filt, [],1);
[r,c]=size(F_filt);
[Cxy1,f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs);
save (['Coherence, misura-C',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Cxy1');

FFT_Kav_gate=FFT_Kav_fft;
for ii=1:length(f)
    if Cxy1(ii)<0.7 & Cxy1(ii)>0.3
        FFT_Kav_gate(ii)=1;
    end
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Confronto PSD calcolata tramite Fft e Psd
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure (109),hold on,

subplot(4,1,[2,3]),hold on
plot (f_fft,20*log10(abs(FFT_Kav_fft).^2),'color',string(colore(i-7,:)),'LineWidth',1)
plot (f(1:length(PSD_Kav_pgram)),20*log10(PSD_Kav_pgram),'.','color',string(colore(i-7,:)),'LineWidth',1)
legend('k(f) tramite Fft (medio sui moduli)','k(f) tramite Psd (medio sui quadrati)')

subplot(4,1,1)
semilogx(f_coer,Cxy1,'color',string(colore(i-7,:)),'LineWidth',2),ylim([0 1])
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

subplot(4,1,[2,3]), hold on,
semilogx (f_fft,20*log10(abs(FFT_Kav_gate).^2),'r.','LineWidth',2)
grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]), %ylim([60 200])

subplot(4,1,4)
semilogx (f_fft,180/pi*angle(FFT_Kav_fft),'color',string(colore(i-7,:)))
grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]);ylim([-180 180])

end
>>>>>>> Stashed changes

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi in intensità PSD
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%calcolo del massimo e del minimo dei plateaux delle PSD in N (quindi
%operando la radice)

bin=round(sqrt(picchi_sel2))+1;
Max_pic=sqrt(max(max(PSD_F))); %calcolo del massimo dei massimi
Min_pic=sqrt(min(max(PSD_F))); %calcolo del minimo dei minimi
delta=(Max_pic-Min_pic)/bin;

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Sostituzione della PSD con la FFT
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<< Updated upstream
PSD_F=fft(F);
PSD_A=fft(A);
f=f2
=======
% PSD_F=2*FFT_F.^2./length (FFT_F(:,1))/fs;
% PSD_A=2*FFT_A.^2./length (FFT_F(:,1))/fs;
% f=f_fft;
>>>>>>> Stashed changes

PSD_F2=PSD_F;
PSD_A2=PSD_A;

PSD_Fsort=[];
PSD_Asort=[];
for i=1:picchi_sel2
    pos=find(sqrt(max(PSD_F2))==sqrt(max(max(PSD_F2))));
    PSD_Fsort=[PSD_Fsort,PSD_F2(:,pos)];
    PSD_Asort=[PSD_Asort,PSD_A2(:,pos)];
    PSD_F2(:,pos)=[];
    PSD_A2(:,pos)=[];
end

%F2=F;A2=A;
figure, hold on
plot(max(PSD_Fsort))
plot(max(PSD_F))
hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi statistica post filtraggio
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

[Y,E] = discretize(sqrt(max(abs(PSD_F))),bin);
values=1:bin;

figure (3), hold on
histfit(sqrt(max(abs(PSD_F))),bin);
hold off
%<<<<<<<<<<<<<<<<<<<<<
% Analisi bin per bin
%<<<<<<<<<<<<<<<<<<<<<

kkk=0;indice=0;
[R,C]=size(PSD_F);
dt=1/fs; time1=1000*(0:dt:r/fs-dt);

for indice = 1:bin
    kkk=kkk+1;
    PSD_Fbin=[];
    PSD_Abin=[];

    FFT_Fbin=[];
    FFT_Abin=[];
    for jj=1:C
        if Y(jj)==indice
            PSD_Fbin=[PSD_Fbin,PSD_F(:,jj)];
            PSD_Abin=[PSD_Abin,PSD_A(:,jj)];
            
            FFT_Fbin=[FFT_Fbin,FFT_F(:,jj)];
            FFT_Abin=[FFT_Abin,FFT_A(:,jj)];
        end
    end
    
    [RR,CC]=size(PSD_Fbin);
    
    if CC>=1

        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo la media degli spettri
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<< Updated upstream
        %PSD
        PSD_Fav =[];
        PSD_Aav =[];
        PSD_Fav = mean(PSD_Fbin, 2);
        PSD_Aav = mean(PSD_Abin, 2);
        PSD_V1av = PSD_Aav./(1i*2*pi*f); %velocità
        PSD_D1av = PSD_V1av./(1i*2*pi*f); % displacement
=======
     
        PSD_Fav_bin = mean(PSD_Fbin,2);
        PSD_Aav_bin = mean(PSD_Abin,2);
        PSD_Vav_bin = PSD_Aav_bin./(1i*2*pi*f).^2; %velocità
        PSD_Dav_bin = PSD_Vav_bin./(1i*2*pi*f).^2; %spostamento
>>>>>>> Stashed changes
        
        FFT_Fav_bin = mean(FFT_Fbin,2);
        FFT_Aav_bin = mean(FFT_Abin,2);
        FFT_Vav_bin = FFT_Aav_bin./(1i*2*pi*f_fft); %velocità
        FFT_Dav_bin = FFT_Vav_bin./(1i*2*pi*f_fft); %spostament
        
        FFT_Kav_bin = FFT_Fav_bin ./ FFT_Dav_bin;
        
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Ciclo for per plottare segnali e spettri (PSD)
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure(101), grid on,
        sgtitle(misura)
        for j=1:C
            if Y(j)==indice
                subplot(2,2,1), hold on
                plot(time1, F(:, j),'color',string(colore(kkk,:))),
                
                subplot(2,2,2), hold on,
                plot(time1, 20*log10(abs(A(:, j))),'color',string(colore(kkk,:)))
                
                subplot(2,2,3), hold on
                semilogx (f, 20*log10(PSD_F(:, j)),'color',string(colore(kkk,:)))
                
                subplot(2,2,4), hold on,
                semilogx (f, 20*log10(PSD_A(:, j)),'color',string(colore(kkk,:)))
            end
        end
        
        subplot(2,2,1), hold on
        xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Force')
        grid on, %xlim([0 10])
        
        subplot(2,2,2), hold on,
        xlabel('Time [ms]'), ylabel('Amplitude 20 log10 ([m/s^2])'), title('Acceleration'),
        grid on, %xlim([0 10])
        
        subplot(2,2,3), hold on
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 N/Hz)'), title('PSD Force')
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on,
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration'),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot degli spettri medi in PSD
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        figure (101), hold on
        
        subplot(2,2,3), hold on,
        plot (f, 20*log10(PSD_Fav_bin), '-.','color',string(colore(kkk,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on,
        plot (f, 20*log10(PSD_Aav_bin), '-.' ,'color',string(colore(kkk,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        saveas (gcf, ['Segnali e spettri-C_',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.fig'])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 20*log10(PSD_Fav_bin);
        fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax=f(fmax(1)+f0);
        %plot sulla PSD della forza
        figure (101)
        subplot (2,2,3),hold on
        xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
        hold off
<<<<<<< Updated upstream

        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo DYNAMIC MASS Mechanical Impedance Dynamic Stiffness
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    %     % Dynamic Mass
    %     DMASS1_av = PSD_Fav./PSD_Aav; %Modulo della Dynamic Mass;
    %     % save (['DMASS1_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwidth),'Hz','.mat'], 'DMASS1_av');
    %     DMASS1_av_ph = angle(FFT_Fav./FFT_Aav); %trovo la fase media usando le FFT;
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

    %     % Mechanical Impedance
    %     MI1_av = PSD_Fav./PSD_V1av; %Modulo dell'Impadenza meccanica
    %     % save (['MI_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwidth),'Hz','.mat'], 'MI1_av');
    %     MI1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_V1av); %trovo la fase media usando le FFT;

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


        %<<<<<<<<<<<<<<<<<<<
        % Dynamic Stiffness
        %<<<<<<<<<<<<<<<<<<<
        Dstiff1_av = PSD_Fav./PSD_D1av; %Modulo della Dynamic Stiffness
        save (['Dstiffness_av, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz - ',num2str(indice),'','.mat'], 'Dstiff1_av');
%         Dstiff1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_D1av); %trovo la fase media usando le FFT;
%         save (['Dstiffness_av_ph, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av_ph');

        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                           
        % Trova picco antirisonanza (frequenza e valore)
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%         y_antiris = max(Dstiff1_av): 
%         pos_y = find (Dstiff1_av == y_antiris);
%         x_antiris = find  
=======
        
        %<<<<<<<<<<<<<<<<<<<
        % Dynamic Stiffness
        %<<<<<<<<<<<<<<<<<<<
>>>>>>> Stashed changes
        
        PSD_Kav_bin = PSD_Fav_bin./PSD_Dav_bin; %Modulo della Dynamic Stiffness
        save (['Dstiffness_av Psd, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz - ',num2str(indice),'','.mat'], 'PSD_Kav_bin');
        %Dstiff1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_D1av); %trovo la fase media usando le FFT;
        %save (['Dstiffness_av_ph, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av_ph');
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % plot singole ripetizioni nel range selezionato
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        figure (round(indice+200))
        hold on
        for iii=1:CC
<<<<<<< Updated upstream
        semilogx (f, 20*log10( PSD_Fbin(:,iii).*((1i*2*pi*f).^2)./PSD_Abin (:,iii) ),'color',string(colore(kkk,:)), 'LineWidth', 1),
        semilogx (f, 20*log10(Dstiff1_av),'k--', 'LineWidth', 1), 
=======
            semilogx (f, 20*log10(abs(PSD_Fbin(:,iii).*((1i*2*pi*f).^4)./PSD_Abin (:,iii) )),'color',string(colore(kkk,:)), 'LineWidth', 1),    
>>>>>>> Stashed changes
        end
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % plot valore medio del bin
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        semilogx (f, 20*log10(PSD_Kav_bin),'k--', 'LineWidth', 2),
        semilogx (f_fft, 20*log10(abs(FFT_Kav_bin.^2)),'r--', 'LineWidth', 2),

        set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
        titolo=['Dynamic Stiffness (Force/Displacement) Amplitude (fron ',num2str(E(indice)),' to ',num2str(E(indice+1)),' N)'];
        title([titolo]),
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
<<<<<<< Updated upstream
        grid on, xlim([ascissamin ascissamax]),ylim([120 220])
=======
        grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
        
>>>>>>> Stashed changes
        saveas (gcf, ['from ',num2str(E(indice)),' to ',num2str(E(indice+1)),' N)',' N Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'.fig'])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot singolo D-Stiffness
        %<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (107),hold on,
        semilogx (f, 20*log10(PSD_Kav_bin),'color',string(colore(kkk,:)), 'LineWidth', 1),
        %semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3),
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        
    end
    
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dstiff totale e settaggio parametri grafico
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (107),  hold on
<<<<<<< Updated upstream
% semilogx (f, 20*log10(Dstiff), 'LineWidth', 3),
grid on, xlim([ascissamin ascissamax]),ylim([120 220])
=======
plot (f_fft, 20*log10(abs(FFT_Kav_fft).^2),'LineWidth',2)
plot (f(1:length(PSD_Kav_pgram)), 20*log10(PSD_Kav_pgram),'.','LineWidth',2)
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
>>>>>>> Stashed changes
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
title(['Dynamic Stiffness (Force/Displacement) Amplitude, Sample: ',campione,'']),

save ('PSD_Dstiff.mat','PSD_Kav_pgram');
saveas (gcf, ['Collezione Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'.fig'])
saveas (gcf, ['Collezione Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'.png'])

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo frequenza di risonanza e K
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
m=1.4;
h=0.005;
s=pi*0.05^2;
lim_sup = find ( f > 1000);
lim_inf = find ( f < 100);
min_k = min(PSD_Kav_pgram(lim_inf(end):lim_sup(1)));
fr = f(find(PSD_Kav_pgram(lim_inf(end):lim_sup(1)) == min_k));
K0=(2*pi*268)^2*m
E=K0*h/s
res=[fr K0 E];
save ('Risultati.mat','res')
<<<<<<< Updated upstream
 figure(107),hold on, plot(f,20*log10(abs(K(1:length(f)))),'r.','LineWidth', 2)
=======
%figure(107),hold on, plot(f,20*log10(abs(K(1:length(f)))),'r.','LineWidth', 2)
>>>>>>> Stashed changes



%%

<<<<<<< Updated upstream
f=f(1:18000);

figure (1),
subplot(3,1,1), 
hold on
semilogx (f, 20*log10(Dstiff1_av_samp1_01_m_1(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp1_01_m_2(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp1_01_m_3(1:18000)),'r', 'LineWidth', 1),

semilogx (f, 20*log10(Dstiff1_av_samp1_02_pl_1(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp1_02_pl_2(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp1_02_pl_3(1:18000)),'b', 'LineWidth', 1),
txt = {'Legenda:', 'Punta: metallo, campione 1, piastra grande, rosso',...
    'Punta: plastica, campione 1, piastra grande, blu'};
 text(120,170,txt,'FontSize',14)
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('log(Frequency) [Hz]','FontSize',10), 
ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]','FontSize', 9), 
title('Punta: plastica, campione 1, piastra grande, biadesivo'), 
grid on
xlim([ascissamin ascissam]), ylim([120 220])
hold off

subplot(3,1,2), hold on
semilogx (f, 20*log10(Dstiff1_av_samp2_01_m_1(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_01_m_2(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_01_m_3(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_1(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_2(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_3(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_4(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_5(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_6(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_7(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_8(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_02_m_9(1:18000)),'r', 'LineWidth', 1),

semilogx (f, 20*log10(Dstiff1_av_samp2_03_pl_1(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_03_pl_2(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_03_pl_3(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_04_pl_1(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_04_pl_2(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_04_pl_3(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_04_pl_4(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_04_pl_5(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp2_04_pl_6(1:18000)),'b', 'LineWidth', 1),
txt = {'Legenda:', 'Punta: metallo, campione 2, piastra grande, rosso',...
    'Punta: plastica, campione 2, piastra grande, blu'};
 text(120,170,txt,'FontSize',14)
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('log(Frequency) [Hz]','FontSize',10), 
ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]','FontSize',9), 
title('Punta: plastica, campione 2, piastra grande, biadesivo'), 
grid on
xlim([100 3000]), ylim([120 220])
hold off

subplot(3,1,3), hold on
semilogx (f, 20*log10(Dstiff1_av_samp3_01_m_1(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_01_m_2(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_01_m_3(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_01_m_4(1:18000)),'r', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_01_m_5(1:18000)),'r', 'LineWidth', 1),

semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_1(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_2(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_3(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_4(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_5(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_6(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_7(1:18000)),'b', 'LineWidth', 1),
semilogx (f, 20*log10(Dstiff1_av_samp3_02_pl_8(1:18000)),'b', 'LineWidth', 1),
txt = {'Legenda:', 'Punta: metallo, campione 3, piastra grande, rosso',...
    'Punta: plastica, campione 3, piastra grande, blu'};
 text(120,170,txt,'FontSize',14)
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('log(Frequency) [Hz]','FontSize',10), 
ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]','FontSize',9), 
title('Punta: plastica, campione 3, piastra grande, biadesivo'), 
grid on
xlim([100 3000]), ylim([120 220])
hold off

saveas (gcf, 'Confronto 2 punte, 3 campioni, piastra grande, biadesivo.fig')

 %%
=======
%%
>>>>>>> Stashed changes

Dstiff1_av=[Dstiff1_av_01, Dstiff1_av_02, Dstiff1_av_03, Dstiff1_av_04, Dstiff1_av_05];
%%
FFT_Kav_pgram =mean(Dstiff1_av, 2);
%%
figure
semilogx (f, 20*log10(Dstiff_media),'b', 'LineWidth', 3),
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('log(Frequency) [Hz]','FontSize',10),
ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]','FontSize',9),
title('Punta: plastica, campione 3, piastra grande, biadesivo'),
grid on
xlim([100 8000]), ylim([120 220])
hold off
