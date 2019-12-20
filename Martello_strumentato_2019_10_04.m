% Martello strumentato filtro intensita
%
set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables

%%
campione={'legno'};
piastra={'piastrina'};
appoggio={'pavimento'};
adesivo={'biadesivo'};
punta={'metallica'};
martellatore=1;
accelerometro=0;
conf = table(campione, piastra, appoggio, adesivo, punta, martellatore, accelerometro)

%%
clc
close all
clear variables

load dati.mat
% load Forza.mat
% load Accelerazione.mat

% conf=[]
%In configuration mettiamo degli indici che ci dicono il campione
%utilizzato, la piastra di carico, la superficie d'appoggio l'adesivo e
%la punta del martello

% campione={'c1'};
% piastra={'pesante1'};
% appoggio={'pavimento'};
% adesivo={'nessuno'};
% punta={'metallica'};
% martellatore=1;
% accelerometro=0;
% conf = table(campione, piastra, appoggio, adesivo, punta, martellatore, accelerometro)

piastre_mass =  [0; 0.006;  0.1967; 0.6274; 0.6249; 1.4293; 2.8871; 15];
piastre_h =     [0; 0.0018; 0.008;  0.008;  0.008;  0.024;  0.0475; 0.16];
piastre_d =     [0; 0.026;  2*sqrt(0.056*0.057/pi); 2*sqrt(0.1*0.1/pi); 2*sqrt(0.1*0.1/pi); 0.1; 0.1;2*sqrt(0.25*0.25/pi)];
piastre = table(piastre_mass,piastre_h,piastre_d);
piastre.Properties.RowNames={'mini','piastrina','quadrata_piccola','quadrata1','quadrata2','pesante1','pesante2','blocco'};
piastre.Properties.VariableNames={'massa','h','d'}

campioni_mass = [0.5531;0.3926; 0.1461  ;0.6128;0.0383;                 0.0705  ;0.1064;0;0;0;0];
campioni_h =    [0.031; 0.027;  0.031   ;0.039; 0.005;                  0.01    ;0.015; 0.045; 0.045; 0.045;0.019];
campioni_d =    [0.1;   0.99;   0.97    ;0.1;   2*sqrt(0.098*0.096/pi); 0.1     ;0.1;   2*sqrt(0.3*0.9/pi);piastre.d(conf.piastra);piastre.d(conf.piastra);piastre.d(conf.piastra)];
campioni = table(campioni_mass,campioni_h,campioni_d);
campioni.Properties.RowNames={'c0','c1','c2','c3','polipropilene','teflon','PVC','slab','viabattelli','viacocchi','legno'};
campioni.Properties.VariableNames={'massa','h','d'}

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
x = piastrina (:,1); % Force [N]
y = piastrina (:,conf.accelerometro+2); % Accelerazione [m/s^2]

% x = reshape(F, [],1);
% y = reshape(A, [],1);

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
L_pre=round((2^15/2)); % Lunghezza della parte prima del picco
L_coda=round(2^15-L_pre-1);     % Lunghezza della coda dei segnali
% Filtraggio doppi colpi
filt_doppi=0;           % Se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi
% Normalizzazione colpi
norm=0;                 % Se norm=1 i colpi vengono normalizzati
% Finestratura
window_F=5;             % Tipo di finestratura da applicare
window_A=5;
wintype = 'hann';
% 0 = nessuna finestratura
% 1 = finestratura quadrata
% 2 = hamming (da verificare)
% 3 = blackmanharris ritardata
% 4 = blackmanharris anticipata
% 5 = Hann Window
% Plotting
ascissamin=1;         % Frequenza minima da plottare nei grafici
ascissamax=5000;       % Frequenza massima da plottare nei grafici
misura = cell2mat(['Campione ',conf.campione,', martellatore ',...
    num2str(conf.martellatore),', punta ',conf.punta,', piastra ',conf.piastra,...
    ',Hann,PSDvsFFT, ',num2str(bandwidth),' Hz']);

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
x =   div_F*x;     % Force [N]
y = g*div_A*(-y);

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
% %F=[]; A=[];
% %load Dati_simulazione
% %A=real(A);
% %Faccio un calcolo di F_filt per ottenere L_win
% [~, L_win] = finestra_forza (F, window_F, fs); 
% %Non applico una vera finestratura.
% [A_filt] = finestra_accelerazione (A, window_A, L_win, fs);
% [F_filt] = finestra_accelerazione (F, window_A, L_win, fs);

i=1 %vuol dire che sto applicando una finestra larga quanto 1/i la lunghezza del segnale.
M_win=r/i; %M_win è l'ordine della finestratura, ossia la sua lunghezza
switch wintype
case 'hann'
win=[zeros((r-M_win)/2,1);hann(M_win);zeros((r-M_win)/2,1)];
case 'rect' %finestratura rettangolare
win=[zeros((r-M_win)/2,1);ones(M_win,1);zeros((r-M_win)/2,1)];
end
E_win=sum(win.^2)/length(win);

%finestratura
A_filt=A.*win;
F_filt=F.*win;
%trasformata (non normalizzata rispetto al numero di punti Fft)
FFT_F=fft(F_filt);
FFT_A=fft(A_filt);
f_fft=0:(r-1);
f_fft=f_fft'/(r-1)*fs;

PSD_F_fft = abs(FFT_F).^2./length (FFT_F(:,1))/fs/E_win; %PSD di F stimata tramite fft
PSD_F_fft(2:end,:)= 2*PSD_F_fft(2:end,:);
PSD_A_fft = abs(FFT_A).^2./length (FFT_A(:,1))/fs/E_win;
PSD_A_fft(2:end,:)= 2*PSD_A_fft(2:end,:);

% %<<<<<<<<<<<<<<<<<<<<<<<<<
% % Abbattimento della coda
% %<<<<<<<<<<<<<<<<<<<<<<<<<
% divid=20;
% for ii=1:c
%     L=L_win(ii); %PER CIASCUNA COLONNA DI f PRENDO LA LUNGHEZZA APPROPRIATA NEL VETTORE LUNGHEZZE L_win
%     F_filt(:,ii)=[F_filt(1:(L-round(L/divid)-1),ii); movavg(F_filt((L-round(L/divid)):end,ii),'linear',round(L/7))];
% end

%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle PSDs
%<<<<<<<<<<<<<<<<<<<<
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
%[Y,E] = discretize(sqrt(max(PSD_F)),bin);
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
plot(f_fft,10*log10(PSD_F_fft(:,1)),'color',string(colore(i,:)))
plot(f,10*log10(PSD_F(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})','PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

subplot (2,1,2),hold on
plot(f_fft,10*log10(PSD_A_fft(:,1)),'color',string(colore(i,:)))
plot(f,10*log10(PSD_A(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})','PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

%<<<<<<<<<<<<<<<<<<<<<
% Calcolo tramite FFT
%<<<<<<<<<<<<<<<<<<<<<

FFT_Fav_fft = mean(FFT_F,2);
FFT_Aav_fft = mean(FFT_A,2);
FFT_Vav_fft = FFT_Aav_fft./(1i*2*pi*f_fft);
FFT_Dav_fft = FFT_Aav_fft./(1i*2*pi*f_fft).^2;
FFT_Kav_fft = FFT_Fav_fft./FFT_Dav_fft;

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% COERENZA usando Forza / Accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

F_filtall = reshape(F_filt, [],1);
A_filtall = reshape(A_filt, [],1);
[r,c]=size(F_filt);
[Cxy1,f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs);
save (cell2mat(['Coerenza_',conf.campione,'_',conf.piastra]), 'Cxy1');

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
sgtitle {'Dynamic Stiffness media K(f)'}

subplot(4,1,[2,3]),hold on
plot (f_fft,10*log10(abs(FFT_Kav_fft).^2),'color',string(colore(i,:)),'LineWidth',1)
plot (f(1:length(PSD_Kav_pgram)),10*log10(PSD_Kav_pgram),'.','color',string(colore(i,:)),'LineWidth',1)
semilogx (f_fft,10*log10(abs(FFT_Kav_gate).^2),'r.','LineWidth',2)
grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]), %ylim([60 200])
legend({'k(f) tramite Fft (medio sui complessi)','k(f) tramite Psd (medio sui quadrati)'},'Location','southeast')

subplot(4,1,1)
semilogx(f_coer,Cxy1,'color',string(colore(i,:)),'LineWidth',2),ylim([0 1])
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

subplot(4,1,4)
semilogx (f_fft,180/pi*angle(FFT_Kav_fft),'color',string(colore(i,:)))
grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]);ylim([-180 180])

subplot(4,1,[2,3]),hold on,ylabel('10\cdotlog_1_0(K(f)) [dB ref. 1 N/m]')
subplot(4,1,[2,3]),hold on,title 'Modulo'
subplot(4,1,[2,3]),ylim([100 200])
subplot(4,1,4),hold on ,ylabel('Angolo [°]'),yticks([-180 -90 0 90 180])
subplot(4,1,1),title 'Coerenza'
subplot(4,1,4),title 'Fase'
xlabel('Frequenza [Hz]')

saveas (gcf, cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'.fig']));
cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'.fig'])

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
% PSD_F=2*FFT_F.^2./length (FFT_F(:,1))/fs;
% PSD_A=2*FFT_A.^2./length (FFT_F(:,1))/fs;
% f=f_fft;

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

m=piastre.massa(conf.piastra); %massa della piastra in uso
h=campioni.h(conf.campione);
s=pi*(campioni.d(conf.campione)/2)^2;
lim_sup = find ( f > 1000);
lim_inf = find ( f < 100);
F_bin=[];
K0_av_bin=[];
E_av_bin=[];
PSD_K_bin=[];
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
        F_bin=[F_bin,(E(kkk)+E(kkk+1))/2];
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo la media degli spettri
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
        PSD_Fav_bin = mean(PSD_Fbin,2);
        PSD_Aav_bin = mean(PSD_Abin,2);
        PSD_Vav_bin = PSD_Aav_bin./(1i*2*pi*f).^2; %velocità
        PSD_Dav_bin = PSD_Vav_bin./(1i*2*pi*f).^2; %spostamento
        
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
                semilogx (f, 10*log10(PSD_F(:, j)),'color',string(colore(kkk,:)))
                
                subplot(2,2,4), hold on,
                semilogx (f, 10*log10(PSD_A(:, j)),'color',string(colore(kkk,:)))
            end
        end
        
        subplot(2,2,1), hold on
        xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Forza')
        grid on, %xlim([0 10])
        
        subplot(2,2,2), hold on,
        xlabel('Time [ms]'), ylabel('Amplitude 20 log10 ([m/s^2])'), title('Accelerazione'),
        grid on, %xlim([0 10])
        
        subplot(2,2,3), hold on
        xlabel('log(Frequency) [Hz]'), ylabel('10 log |PSD| (dB ref 1 N/Hz)'), title('PSD Forza')
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on,
        xlabel('log(Frequency) [Hz]'), ylabel('10 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Accelerazione'),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot degli spettri medi in PSD
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        figure (101), hold on
        
        subplot(2,2,3), hold on,
        plot (f, 10*log10(PSD_Fav_bin), '-.','color',string(colore(kkk,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on,
        plot (f, 10*log10(PSD_Aav_bin), '-.' ,'color',string(colore(kkk,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        saveas (gcf, ['Segnali e spettri.fig'])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 10*log10(PSD_Fav_bin);
        fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax=f(fmax(1)+f0);
        %plot sulla PSD della forza
        figure (101)
        subplot (2,2,3),hold on
        xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
        hold off
        
        %<<<<<<<<<<<<<<<<<<<
        % Dynamic Stiffness
        %<<<<<<<<<<<<<<<<<<<
        
        PSD_Kav_bin = PSD_Fav_bin./PSD_Dav_bin; %Modulo della Dynamic Stiffness
        save (['Dstiffness_av_bin.mat'], 'PSD_Kav_bin');
        %Dstiff1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_D1av); %trovo la fase media usando le FFT;
        %save (['Dstiffness_av_ph, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av_ph');
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % plot singole ripetizioni nel range selezionato
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        figure (round(indice+200))
        hold on
        for iii=1:CC
            semilogx (f, 10*log10(abs(PSD_Fbin(:,iii).*((1i*2*pi*f).^4)./PSD_Abin (:,iii) )),...
                'color',string(colore(kkk,:)), 'LineWidth', 1),    
        end
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % plot valore medio del bin
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        semilogx (f,    10*log10(PSD_Kav_bin),'k--', 'LineWidth', 2),
        semilogx (f_fft,10*log10(abs(FFT_Kav_bin.^2)),'r--', 'LineWidth', 2),

        set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
        
        titolo=['Dynamic Stiffness (Force/Displacement) Amplitude (fron ',num2str(E(indice)),...
            ' to ',num2str(E(indice+1)),' N)'];
        
        title([titolo]),
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
        
        saveas (gcf, ['from ',num2str(E(indice)),' to ',num2str(E(indice+1)),' N)',' N Dstiff.fig'])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot singolo D-Stiffness
        %<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (107),hold on,
        semilogx (f, 10*log10(PSD_Kav_bin),'color',string(colore(kkk,:)), 'LineWidth', 1),
        %semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3),
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        

        min_k = min(PSD_Kav_bin(lim_inf(end):lim_sup(1)));
        fr = f(lim_inf(end)+find(PSD_Kav_bin(lim_inf(end):lim_sup(1)) == min_k));
        K0_av_bin=[K0_av_bin,(2*pi*fr)^2*m];
        E_av_bin=[E_av_bin,K0_av_bin*h/s];
        PSD_K_bin = [PSD_K_bin,K0_av_bin];

    end
    
end

save (cell2mat(['Dstiffness_bin_',conf.campione,'_',conf.piastra,'_Fft.mat']),'PSD_K_bin');


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dstiff totale e settaggio parametri grafico
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (107),  hold on
plot (f_fft, 10*log10(abs(FFT_Kav_fft).^2),'LineWidth',2)
plot (f(1:length(PSD_Kav_pgram)), 10*log10(PSD_Kav_pgram),'.','LineWidth',2)
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('Frequenza [Hz]'), ylabel('Psd Rigidezza Dinamica [dB ref 1 N/m]'),
title(['Pds della rigidezza dinamica [N/mHz]. Campione: ',cell2mat([conf.campione,...
    ' + ',conf.piastra,'. Adesivo: ',conf.adesivo,'.'])]),

save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Psd.mat']),'PSD_Kav_pgram');
save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Fft.mat']),'FFT_Kav_fft');

saveas (gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.fig']))
saveas (gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.png']))

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo frequenza di risonanza e K
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

m=piastre.massa(conf.piastra); %massa della piastra in uso
h=campioni.h(conf.campione);
s=pi*(campioni.d(conf.campione)/2)^2;
lim_sup = find ( f > 1000);
lim_inf = find ( f < 100);
min_k = min(PSD_Kav_pgram(lim_inf(end):lim_sup(1)));
fr = f(lim_inf(end)+find(PSD_Kav_pgram(lim_inf(end):lim_sup(1)) == min_k));
K0_av=(2*pi*fr)^2*m
E_av=K0_av*h/s

%figure(107),hold on, plot(f,20*log10(abs(K(1:length(f)))),'r.','LineWidth', 2)


%<<<<<<<<<<<<<
% Stampo K(F)
%<<<<<<<<<<<<<

figure(120)
plot(F_bin,K0_av_bin,'+')
P = polyfit(F_bin,K0_av_bin,1);
hold on
plot(F_bin,F_bin*P(1)+P(2),'-')
grid on
saveas (gcf, cell2mat(['Fit_Dstiffness_',conf.campione,'_',conf.piastra,'.fig']))

K0_fit=P(2)
E_fit=K0_fit*h/s
risultati = table(fr,K0_av,E_av,K0_fit,E_fit)
save (cell2mat(['Risultati_',conf.campione,'_',conf.piastra,'.mat']),'risultati')

risultati_bin = table(F_bin',K0_av_bin');
risultati_bin.Properties.VariableNames={'F_bin','K0_av_bin'}
save (cell2mat(['Risultati_bin_',conf.campione,'_',conf.piastra,'.mat']),'risultati_bin')

%%

%%

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
