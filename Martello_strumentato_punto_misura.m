% Martello strumentato
% Versione per elzborare una singola registrazione
%
set (0,'DefaultFigureWindowStyle','docked')
clc
close all

%%
clear variables

%%
conf=[]

%In configuration mettiamo degli indici che ci dicono il campione
%utilizzato, la piastra di carico, la superficie d'appoggio l'adesivo e
%la punta del martello

campione={'polipropilene'};
piastra={'pesante1'};
appoggio={'nessuno'};
adesivo={'nessuno'};
punta={'metallica'};
martellatore=1;
accelerometro=0;
conf = table(campione, piastra, appoggio, adesivo, punta, martellatore, accelerometro)

%%
clc
close all
clear variables

load dati1.mat

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf,piastre);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
x1 = data (:,1); % Force [N]
y1 = data (:,conf.accelerometro+2); % Accelerazione [m/s^2]
clear data;

load dati2.mat
x2 = data (:,1); % Force [N]
y2 = data (:,conf.accelerometro+2); % Accelerazione [m/s^2]
clear data;

load dati3.mat
x3 = data (:,1); % Force [N]
y3 = data (:,conf.accelerometro+2); % Accelerazione [m/s^2]
clear data;

%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri di controllo 
%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri fisici
[g,div_F,div_A,fs] = parametri_fisici();

% Parametri di ricerca
[soglia,delay,inizio,fine1] = parametri_ricerca_picchi(fs,x1);
[soglia,delay,inizio,fine2] = parametri_ricerca_picchi(fs,x2);
[soglia,delay,inizio,fine3] = parametri_ricerca_picchi(fs,x3);
% Parametri di creazione dei campioni
[Lsample,L_pre,L_coda] = parametri_creamatrice();

% Parametri di filtro
bandwidth=0;
filt_doppi=0;           % Se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi

% Finestratura
wintype = 'hann';
% Determina il tipo di finestratura da utilizzare per i segnali di
% accelerazione:
% hann = utilizza la finestra di hann;
% rect = utilizza una finestratura rettangolare;
% none = non applica nessuna finestratura.

% Plotting
ascissamin=80;          % Frequenza minima da plottare nei grafici
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
x1 =   div_F*x1;     % Force [N]
y1 = g*div_A*(-y1);
x2 =   div_F*x2;     % Force [N]
y2 = g*div_A*(-y2);
x3 =   div_F*x3;     % Force [N]
y3 = g*div_A*(-y3);

% Plot dei segnali
L = length(x1);
dt=1/fs; time=[0:dt:L/fs-dt];
figure(1)
subplot(2,1,1), hold on, plot (time, x1)
subplot(2,1,2), hold on, plot (time, y1),
hold off

L = length(x2);
dt=1/fs; time=[0:dt:L/fs-dt];
figure(2)
subplot(2,1,1), hold on, plot (time, x2)
subplot(2,1,2), hold on, plot (time, y2),
hold off

L = length(x3);
dt=1/fs; time=[0:dt:L/fs-dt];
figure(3)
subplot(2,1,1), hold on, plot (time, x3)
subplot(2,1,2), hold on, plot (time, y3),
hold off

%<<<<<<<<<<<<<<<<<<<<
% Ricerca dei PICCHI
%<<<<<<<<<<<<<<<<<<<<
% inizio=29.3*fs;
[picchi1,n_picchi1] = trovacolpi(x1, soglia, delay, inizio, fine1);
n_picchi1
% inizio=0.1*fs
[picchi2,n_picchi2] = trovacolpi(x2, soglia, delay, inizio, fine2);
n_picchi2

[picchi3,n_picchi3] = trovacolpi(x3, soglia, delay, inizio, fine3);
n_picchi3

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[F1, pos] = creamatriceforza_noavg (x1, picchi1,n_picchi1, L_pre, L_coda, filt_doppi, fs);
[A1] = creamatriceaccelerazione (y1, pos, L_pre, L_coda, fs);
picchi_sel1=length(pos)

pos=[];
[F2, pos] = creamatriceforza_noavg (x2, picchi2,n_picchi3, L_pre, L_coda, filt_doppi, fs);
[A2] = creamatriceaccelerazione (y2, pos, L_pre, L_coda, fs);
picchi_sel2=length(pos)

pos=[];
[F3, pos] = creamatriceforza_noavg (x3, picchi3,n_picchi3, L_pre, L_coda, filt_doppi, fs);
[A3] = creamatriceaccelerazione (y3, pos, L_pre, L_coda, fs);
picchi_sel3=length(pos)

% Accodo tutti i segnali di forza e accelerazione in un'unica matrice
F=[F1 F2 F3];
A=[A1 A2 A3];

% Matrice degli indici
% Se l'elemento i-esimo di I vale x, la colonna i-esima di F appartiene
% alla matrice forza Fx.
I=[ones(size(F1(1,:))) 2*ones(size(F2(1,:))) 3*ones(size(F3(1,:)))];
I0=I;
[r,c]=size(F);

%<<<<<<<<<<<<<<
% Finestratura
%<<<<<<<<<<<<<<
% Generazione finestra Accelerazione
M_win=L_pre; %M_win è l'ordine della finestratura, ossia la sua lunghezza
switch wintype
    case 'hann'
        curve=hann(M_win);
        win_A=[curve(1:end/2);ones((r-M_win-1),1);curve(end/2:end)];
    case 'rect' %finestratura rettangolare
        curve=ones(M_win);
        win_A=[curve(1:end/2);ones((r-2*M_win-1),1);curve(end/2:end)];
    case 'none'
        win_A=ones(r,1);
end
E_win=sum(win_A.^2)/length(win_A);
dt=1/fs; time1=1000*(0:dt:r/fs-dt);

% Generazione finestra Forza
L_win_F=L_pre+round(1*fs/1000)+M_win/2;
win_F=[ones(L_win_F-(M_win/2),1);curve(end/2:end-1);zeros(Lsample-L_win_F,1)];
win_F(1,1)=0;

% Plot delle finestre
figure(101),hold on
subplot(2,2,2)
plot(time1,win_A*20*log10(max(max(A))));
subplot(2,2,1);
plot(time1,win_F*max(max(F1)));
hold off

% Finestratura
F_filt=F.*win_F;
A_filt=A.*win_A;

%<<<<<<<<<<<<<<<<<<<<<
% Analisi smorzamento
%<<<<<<<<<<<<<<<<<<<<<
A_reverse=cumsum(flip(A.^2));
figure (14), hold on,
plot(10*log10(movmean(A.^2,20)))
figure(15), hold on, plot(10*log10(flip(A_reverse)))

%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle PSDs
%<<<<<<<<<<<<<<<<<<<<
% Trasformata (non normalizzata rispetto al numero di punti Fft)
FFT_F=fft(F_filt);
FFT_A=fft(A_filt);

% Frequenza della FFT_F
f_fft=0:(r-1);
f_fft=f_fft'/(r-1)*fs;

% Calcolo della PSD tramite FFT
PSD_F_fft = abs(FFT_F).^2./length(FFT_F(:,1))/fs/1; %PSD di F stimata tramite fft
PSD_F_fft(2:end,:)= 2*PSD_F_fft(2:end,:);
PSD_A_fft = abs(FFT_A).^2./length (FFT_A(:,1))/fs/E_win;
PSD_A_fft(2:end,:)= 2*PSD_A_fft(2:end,:);

% Calcolo della PSD tramite Periodogram
win_1=ones(size(win_F)); % finestra unitaria per non far calcolare la normalizzazione a periodogram
PSD_F = []; PSD_A = [];
% Periodogram non gestisce più di 40 colonne contemporaneamente quindi
% eseguo le operazioni ciclicamente
for i=1:3
    [PSD_Ftemp, f]= periodogram(F_filt(:,I==i), win_1, r, fs); %PSD Forza [N^2]
    [PSD_Atemp, f]= periodogram(A(:,I==i), win_A, r, fs); %PSD Accelerazione [g^2]
    PSD_F=[PSD_F PSD_Ftemp];
    PSD_A=[PSD_A PSD_Atemp];
end
clear PSD_Ftemp
clear PSD_Atemp
save ('f.mat', 'f');

%<<<<<<<<<<<<<<<<<<<<<<
% Filtraggio per Banda
%<<<<<<<<<<<<<<<<<<<<<<
[R,C]=size(PSD_F);
tagli=[];
scarti=0;
f0=find(f>ascissamin,1);
for jj=1:C
    fmax=find(PSD_F(f0:end, jj)<((PSD_F(f0, jj)/10)),1);
    fmax=f(fmax+f0);
    if  fmax<bandwidth
        tagli=[tagli; jj];
    end
end
scarti=length(tagli);
picchi_sel2 = picchi_sel1 - scarti

% Eliminazione degli scarti
% PSD_F(:,tagli)=[]; % Cancello le forze
% PSD_A(:,tagli)=[]; % Cancello le accelerazioni
% PSD_F_fft(:,tagli)=[]; % Cancello le forze fft
% PSD_A_fft(:,tagli)=[]; % Cancello le accelerazioni fft
I(:,tagli)=0; % Cancello gli indici corrispondenti

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi statistica pre filtraggio
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
bin=[];
Y=[];
E=[];
% bin=round(sqrt(picchi_sel2))+1;
% %[Y,E] = discretize(sqrt(max(PSD_F)),bin);
% values=1:bin;
% figure (3), hold on
% histfit(sqrt(max(PSD_F)),bin);

%<<<<<<<<<<<<<<<<<<<<<<<
% Filtraggio Intensita'
%<<<<<<<<<<<<<<<<<<<<<<<
[~,c]=size(PSD_F);
tagli=[];
scarti=0;
for jj=1:(c)
    if (max(PSD_F(:,jj))< 0.1*max(max(PSD_F)))
        tagli=[tagli; jj];
    end
end
picchi_sel2 = picchi_sel2 - scarti
% PSD_F(:,tagli)=[];
% PSD_A(:,tagli)=[];
% PSD_F_fft(:,tagli)=[];
% PSD_A_fft(:,tagli)=[];
I(:,tagli)=0;

%<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo Dstiff totale
%<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo tramite periodogram
PSD_Fav_pgram = mean(PSD_F(:,I~=0), 2);
PSD_Aav_pgram = mean(PSD_A(:,I~=0), 2);
PSD_Vav_pgram = PSD_Aav_pgram./(1i*2*pi*f).^2; % velocitÃ 
PSD_Dav_pgram = PSD_Vav_pgram./(1i*2*pi*f).^2; % displacement
PSD_Kav_pgram = PSD_Fav_pgram./PSD_Dav_pgram; % dynamic stiffness calcolata tramite periodogram

% Calcolo delle FFT (serve per la fase)
% FFT_Fav_fft = mean(FFT_F,2);
% FFT_Aav_fft = mean(FFT_A,2);
% FFT_Vav_fft = FFT_Aav_fft./(1i*2*pi*f_fft);
% FFT_Dav_fft = FFT_Aav_fft./(1i*2*pi*f_fft).^2;
% FFT_Kav_fft = FFT_Fav_fft./FFT_Dav_fft;

% Plot Forza media vs Accelerazione media e confronto tra periodogram e fft
% i=1;
% figure(110), hold on
%
% subplot (2,1,1),hold on,
% plot(f_fft,10*log10(mean(PSD_F_fft,2)),'color',string(colore(i,:)))
% plot(f,10*log10(PSD_Fav_pgram(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
% grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
% legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})',...
%     'PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')
%
% subplot (2,1,2),hold on
% plot(f_fft,10*log10(mean(PSD_A_fft,2)),'color',string(colore(i,:)))
% plot(f,10*log10(PSD_Aav_pgram(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
% grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
% legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})',...
%     'PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')
%
% saveas(gcf,'Forza_vs_Accelerazione.fig');

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% COERENZA usando Forza / Accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo le tre coerenze. forse è il caso di spostare tutto nel ciclo
% principale
% for i=1:3
%     % Accodo tutti i segnali della stessa misura
%     F_filtall = reshape(F_filt(:,I==i), [],1);
%     A_filtall = reshape(A_filt(:,I==i), [],1);
%     [r,c]=size(F_filt);
%     % Calcolo coerenza
%     [Cxy(:,i),f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs);
% end
% save (cell2mat(['Coerenza_',conf.campione,'_',conf.piastra]), 'Cxy1');

% PSD_Kav_gate=PSD_Kav_pgram;
% f_Kav_gate=f;
% tagli=[];
% for ii=1:length(f)
%     if Cxy(ii)<0.7 & Cxy(ii)>0.3
%         tagli=[tagli,ii];
%     end
% end
% PSD_Kav_gate(tagli)=[];
% f_Kav_gate(tagli)=[];

% Acc=PSD_Aav_pgram./PSD_Fav_pgram;
% save(('Acc'), 'Acc');
% save(('Forza'), 'PSD_Fav_pgram');
% save(('Accelerazione'), 'PSD_Aav_pgram');
%
% figure(111), hold on
%
% subplot (2,1,1),hold on,
% plot(f_coer,Cxy,'LineWidth',2)
% grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
% legend('Coerenza')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')
%
% subplot (2,1,2),hold on
% plot(f,Acc(1:length(f)),'LineWidth',2)
% grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
% legend('Accelerazione su Forza')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo della f massima media
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% PSD_Fav_dB = 10*log10(PSD_Fav_pgram);
% fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
% fmax=f(fmax(1)+f0);
% %plot sulla PSD della forza
% figure (111)
% subplot (2,1,2),hold on
% xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
% hold off
% saveas(gcf,'Acc.fig');
%
% fr=f(find (Acc==max(Acc(find (f>ascissamin,1):find (f>fmax,1)))));
%
% m=piastre.massa(conf.piastra); %massa della piastra in uso
% h=campioni.h(conf.campione);
% s=pi*(campioni.d(conf.campione)/2)^2;
%
% K0_av=(2*pi*fr)^2*m
% St=K0_av/s
% E_av=K0_av*h/s
%
% Risultati_Acc = table(fr,K0_av,St,E_av);
% Risultati_Acc.Properties.VariableNames={'f_r','K_0','S_t','E'}
% save (['Risultati_Acc.mat'], 'Risultati_Acc');

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Confronto PSD calcolata tramite Fft e Psd
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% figure (109),hold on,
% sgtitle {'Dynamic Stiffness media K(f)'}
%
% subplot(4,1,[2,3]),hold on
% plot (f_fft,10*log10(abs(FFT_Kav_fft).^2),'color',string(colore(i,:)),'LineWidth',1)
% plot (f(1:length(PSD_Kav_pgram)),10*log10(PSD_Kav_pgram),'.','color',string(colore(i,:)),'LineWidth',1)
% semilogx (f_Kav_gate,10*log10(abs(PSD_Kav_gate).^2),'r.','LineWidth',2)
% grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]), %ylim([60 200])
% legend({'k(f) tramite Fft (medio sui complessi)','k(f) tramite Psd (medio sui quadrati)'},'Location','southeast')
%
% subplot(4,1,1)
% semilogx(f_coer,Cxy1,'color',string(colore(i,:)),'LineWidth',2),ylim([0 1])
% grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
%
% subplot(4,1,4)
% semilogx (f_fft,180/pi*angle(FFT_Kav_fft),'color',string(colore(i,:)))
% grid on, set(gca, 'XScale', 'log'),xlim([ascissamin ascissamax]);ylim([-180 180])
%
% subplot(4,1,[2,3]),hold on,ylabel('10\cdotlog_1_0(K(f)) [dB ref. 1 N/m]')
% subplot(4,1,[2,3]),hold on,title 'Modulo'
% subplot(4,1,[2,3]),ylim([100 200])
% subplot(4,1,4),hold on ,ylabel('Angolo [Â°]'),yticks([-180 -90 0 90 180])
% subplot(4,1,1),title 'Coerenza'
% subplot(4,1,4),title 'Fase'
% xlabel('Frequenza [Hz]')
%
% saveas (gcf, cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'.fig']));
% cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'.fig'])

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi in intensitÃ  PSD
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% calcolo del massimo e del minimo dei plateaux delle PSD in N (quindi
% operando la radice)
% bin=round(sqrt(picchi_sel2))+1;
% Max_pic=sqrt(max(max(PSD_F))); %calcolo del massimo dei massimi
% Min_pic=sqrt(min(max(PSD_F))); %calcolo del minimo dei minimi
% delta=(Max_pic-Min_pic)/bin;

%<<<<<<<<<<<<<<<<<<<<<<<
% Plot del massimi di F
%>>>>>>>>>>>>>>>>>>>>>>>
% PSD_F2=PSD_F;
% PSD_A2=PSD_A;
% % Ricerca dei massimi e ordinamento
% PSD_Fsort=[];
% PSD_Asort=[];
% for i=1:picchi_sel2
%     pos=find(sqrt(max(PSD_F2))==sqrt(max(max(PSD_F2))));
%     PSD_Fsort=[PSD_Fsort,PSD_F2(:,pos)];
%     PSD_Asort=[PSD_Asort,PSD_A2(:,pos)];
%     PSD_F2(:,pos)=[];
%     PSD_A2(:,pos)=[];
% end
% % plot dei massimi
% figure (5), hold on
% ylimsup =  floor(max(max(PSD_F)));
% ylim=([0 1]);
% plot(max(PSD_Fsort));
% plot(max(PSD_F));
% hold off


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi statistica post filtraggio
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% [Y,E] = discretize(sqrt(max(abs(PSD_F))),bin);
% values=1:bin;
%
% figure (3), hold on
% histfit(sqrt(max(abs(PSD_F))),bin);
% hold off

%<<<<<<<<<<<<<<<<<<<<<
% Analisi bin per bin
%<<<<<<<<<<<<<<<<<<<<<
m=piastre.massa(conf.piastra); %massa della piastra in uso
h=campioni.h(conf.campione);
s=pi*(campioni.d(conf.campione)/2)^2;
lim_sup = find ( f > 1000);
lim_inf = find ( f < 100);
K0_av_bin=[];
E_av_bin=[];
PSD_K_bin=[];
[R,C]=size(PSD_F);
F_max=[];
% Parametri di ingresso:
% PSD_F, PSD_A, F_filt, A_filt,I
for mis = 1:3
    
    % Seleziono i campioni appartenenti al BIN
    PSD_F_misura = PSD_F(:,I==mis);
    PSD_A_misura = PSD_A(:,I==mis);
    
    % Calcolo dimensioni del BIN
    [RR,CC]=size(PSD_F_misura);
    
    if CC>=1 % se il BIN non è vuoto

        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo Forza media di ciascuna registrazione
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        F_max_av (mis)=mean(max(F(:,I==mis)));
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % COERENZA usando Forza / Accelerazione
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo le tre coerenze. forse è il caso di spostare tutto nel ciclo
        % principale
        F_filtall = reshape(F_filt(:,I==mis), [],1);
        A_filtall = reshape(A_filt(:,I==mis), [],1);
        [r,c]=size(F_filt);
        
        % Calcolo coerenza
        [Cxy(:,mis),f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs);
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo la media degli spettri
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_misura(:,mis) = mean(PSD_F_misura,2);
        PSD_Aav_misura(:,mis) = mean(PSD_A_misura,2);
        PSD_Vav_misura(:,mis) = PSD_Aav_misura(:,mis)./(1i*2*pi*f).^2; %velocitÃ 
        PSD_Dav_misura(:,mis) = PSD_Vav_misura(:,mis)./(1i*2*pi*f).^2; %spostamento
        % Dynamic Stiffness
        PSD_Kav_misura(:,mis) = PSD_Fav_misura(:,mis)./PSD_Dav_misura(:,mis); %Modulo della Dynamic Stiffness
        
        
        
        FFT_F_misura(:,mis) = mean(FFT_F(:,I==mis),2);
        FFT_A_misura(:,mis) = mean(FFT_A(:,I==mis),2);
        FFT_V_misura(:,mis) = FFT_A_misura(:,mis)./(1i*2*pi*f_fft);
        FFT_D_misura(:,mis) = FFT_A_misura(:,mis)./(1i*2*pi*f_fft).^2;
        FFT_K_misura(:,mis) = FFT_F_misura(:,mis)./FFT_D_misura(:,mis);
        
        
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo deviazione standard
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        devst_F_misura(:,mis) = std (abs(PSD_F_misura),0,2);
        devst_A_misura(:,mis) = std (abs(PSD_A_misura),0,2);
        
        dF = devst_F_misura(:,mis);
        dA  = devst_A_misura(:,mis);
        b = abs((1i*2*pi*f).^4);
        FF = abs(PSD_Fav_misura(:,mis));
        AA = abs(PSD_Aav_misura(:,mis));
        C = abs(Cxy(1:1+end/2,mis));
        KK = PSD_Kav_misura(:,mis);
        % K = F*(1i*2*pi).^4 /A
        
        dK = ( ((b./AA).*dF).^2 + ((-b.*FF./AA.^2).*dA).^2 + C.*(b./AA).*(-b.*FF./AA.^2).*dF.*dA  ).^(1/2);
        devst_K_misura(:,mis) = dK;
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 10*log10(PSD_Fav_misura(:,mis));
        fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax=f(fmax(1)+f0);
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot di segnali e spettri (PSD)
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure(101), grid on,
        sgtitle(misura)
        for i=1:c
            if I(i)==mis
                subplot(2,2,1), hold on
                plot(time1, F(:,i),'color',string(colore(mis,:))),
                
                subplot(2,2,2), hold on,
                plot(time1, 20*log10(abs(A(:, i))),'color',string(colore(mis,:)))
                
                subplot(2,2,3), hold on
                semilogx (f, 10*log10(PSD_F(:, i)),'--','color',string(colore(mis,:)))
                
                subplot(2,2,4), hold on,
                semilogx (f, 10*log10(PSD_A(:, i)),'--','color',string(colore(mis,:)))
            end
        end
        
        % Plot spettri medi
        subplot(2,2,3), hold on,
        plot (f, 10*log10(PSD_Fav_misura(:,mis)),'color',string(colore(mis,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        subplot(2,2,4), hold on,
        plot (f, 10*log10(PSD_Aav_misura(:,mis)),'color',string(colore(mis,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        %plot frequenza massima sulla PSD della forza
        subplot (2,2,3), hold on
        xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
        hold off
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot grafico riassuntivo con dstiff medi di tutti le misure
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (107),hold on,
        
        subplot(4,1,1), hold on
        % plot coerenza
        semilogx (f, C,'color',string(colore(mis,:)), 'LineWidth', 2),
        
        subplot(4,1,[2,3]), hold on
        % Plot rigidezza dinamica
        semilogx (f, 10*log10(PSD_Kav_misura(:,mis)),'color',string(colore(mis,:)), 'LineWidth', 2),
        % Plot deviazione standard
        semilogx (f, 10*log10(PSD_Kav_misura(:,mis)-dK),'-.','color',string(colore(mis,:)), 'LineWidth', 1),
        semilogx (f, 10*log10(PSD_Kav_misura(:,mis)+dK),'-.','color',string(colore(mis,:)), 'LineWidth', 1),
        
        % Plot limite in frequenza
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(mis,:)));
        
        % Plot fase
        subplot(4,1,4), hold on
        semilogx (f_fft,180/pi*angle(FFT_K_misura(:,mis)),'color',string(colore(mis,:)),'LineWidth', 2)
        
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot Accelerazione e forza
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<        
        figure(110), hold on
        subplot (2,1,1),hold on,
        plot(f,10*log10(FF),'color',string(colore(mis,:)),'LineWidth',2)
        plot(f,10*log10(FF+dF),'--','color',string(colore(mis,:)),'LineWidth',1)
        plot(f,10*log10(FF-dF),'--','color',string(colore(mis,:)),'LineWidth',1)
        
        subplot (2,1,2),hold on,
        plot(f,10*log10(AA),'color',string(colore(mis,:)),'LineWidth',2)
        plot(f,10*log10(AA+dA),'--','color',string(colore(mis,:)),'LineWidth',1)
        plot(f,10*log10(AA-dA),'--','color',string(colore(mis,:)),'LineWidth',1)
        
        % plot(f,10*log10(PSD_Fav_pgram(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
        % grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
        % legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})',...
        %     'PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')
        %
        % subplot (2,1,2),hold on
        % plot(f_fft,10*log10(mean(PSD_A_fft,2)),'color',string(colore(i,:)))
        % plot(f,10*log10(PSD_Aav_pgram(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
        % grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
        % legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})',...
        %     'PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')
        %
        % saveas(gcf,'Forza_vs_Accelerazione.fig');
        
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo K0 della misura
        %<<<<<<<<<<<<<<<<<<<<<<<<<
        min_k = min(PSD_Kav_misura(lim_inf(end):lim_sup(1),mis));
        fr = f(lim_inf(end)+find(PSD_Kav_misura(lim_inf(end):lim_sup(1),mis) == min_k));
        K0_av_bin=[K0_av_bin,(2*pi*fr)^2*m];
        E_av_bin=[E_av_bin,K0_av_bin*h/s];
        PSD_K_bin = [PSD_K_bin,K0_av_bin];
        
    end
    
end
% Salvataggio coerenza
save (cell2mat(['Coerenza_',conf.campione,'_',conf.piastra]), 'Cxy');

% Salvataggio rigidezza dinamica
save (['Dstiffness_av_misura.mat'], 'PSD_Kav_misura');

% Figura Segnali e Spettri
% Settaggio parametri grafico
figure (101), hold on
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
% Salvataggio
figure (101),saveas (gcf, ['Segnali e spettri.fig'])

% Salvataggio delle k0 misurate bin per bin
save (cell2mat(['Dstiffness_bin_',conf.campione,'_',conf.piastra,'_Fft.mat']),'PSD_K_bin');

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dstiff totale e settaggio parametri grafico
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (107)
sgtitle(['PSD della rigidezza dinamica [N/mHz]. Campione: ',cell2mat([conf.campione,...
    ' + ',conf.piastra,'. Adesivo: ',conf.adesivo,'.'])])

subplot(4,1,1), hold on
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),ylim([0 1])
xlabel('Frequenza [Hz]'), ylabel('Coerenza'),
legend(['misura 1 (',num2str(round(F_max_av (1))),' N)'],...
    ['misura 2 (',num2str(round(F_max_av (2))),' N)'],...
    ['misura 3 (',num2str(round(F_max_av (3))),' N)'])

% Plot della massa della piastra
subplot(4,1,[2,3]), hold on
plot(f,10*log10(m*(2*pi*f).^4),'--k');
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Rigidezza Dinamica [dB @ 1 N/m]'),

subplot(4,1,4), hold on
set(gca, 'XScale', 'log'), yticks([-180 -90 0 90 180])
grid on, xlim([ascissamin ascissamax]),ylim([-180 180])
xlabel('Frequenza [Hz]'), ylabel('Fase [ang]'),
legend(['misura 1 (',num2str(round(F_max_av (1))),' N)'],...
    ['misura 2 (',num2str(round(F_max_av (2))),' N)'],...
    ['misura 3 (',num2str(round(F_max_av (3))),' N)'])

% Salvataggio
saveas (gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.fig']))
saveas (gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.png']))



%<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Accelerazione e forza
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure(110), hold on

subplot (2,1,1),hold on,
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Forza [dB @ 1 N]'),

subplot (2,1,2),hold on,
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Accelerazione [dB @ 1 m/s^2]'),
saveas (gcf, cell2mat(['Forza_vs_Accelerazione_',conf.campione,'_',conf.piastra,'.fig']))


% Salvataggio delle dstiffness medie FFT e PSD
% save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Psd.mat']),'PSD_Kav_pgram');
% save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Fft.mat']),'FFT_Kav_fft');


return

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

%%


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo deviazione standard
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
devst_F_misura(:,mis)= std (abs(PSD_F_misura),0,2);
devst_A_misura(:,mis)= std (abs(PSD_A_misura),0,2);

dF = devst_F_misura(:,mis);
dA  = devst_A_misura(:,mis);
b = abs((1i*2*pi*f).^4);
F = abs(PSD_Fav_misura(:,mis));
A = abs(PSD_Aav_misura(:,mis));
C = abs(Cxy(1:1+end/2,mis));
K = PSD_Kav_misura(:,mis);
% K = F*(1i*2*pi).^4 /A

dK= ( ((b./A).*dF).^2 + ((-b.*F./A.^2).*dA).^2 + C.*(b./A).*(-b.*F./A.^2).*dF.*dA  ).^(1/2);

% k= abs (F b /A) 
figure, hold on
plot (f, 10*log10(F))
plot (f, 10*log10(dF))
plot (f, 10*log10(A))
plot (f, 10*log10(dA))
plot (f, 10*log10(K))
plot (f, 10*log10(dK))

plot (f, 10*log10(abs(K)+dK),'--')
plot (f, 10*log10(abs(K)-dK),'--')

legend('F','dF','A','dA','K','dK','K+dK','K-dK')
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
%%
% TEST CALCOLO PSD

win_1=ones(size(win_F));
[PSD_F, f]= periodogram(F, win_F, r, fs); %PSD Forza [N^2]
[PSD_F1, f]= periodogram(F, win_1, r, fs); %PSD Forza [N^2]
[PSD_F2, f]= periodogram(F_filt, win_F, r, fs); %PSD Forza [N^2]
[PSD_F3, f]= periodogram(F_filt, win_1, r, fs); %PSD Forza [N^2]

figure(56),hold on
plot(f,20*log10(PSD_F(:,1)));plot(f,20*log10(PSD_F1(:,1)))
plot(f,20*log10(PSD_F2(:,1))),plot(f,20*log10(PSD_F3(:,1)))
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])

legend('F win_f','F win_1','F_{filt} win_f','F_{filt} win_1')

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
