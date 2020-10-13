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

campione={'c2'};
piastra={'pesante1'};
appoggio={'nessuno'};
adesivo={'gesso'};
punta={'gomma'};
martellatore=1;
accelerometro=0;
conf = table(campione, piastra, appoggio, adesivo, punta, martellatore, accelerometro)

%%
clc
close all
clear variables

load dati.mat

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf,piastre);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
x = data (:,1); % Force [N]
y = data (:,conf.accelerometro+2); % Accelerazione [m/s^2]

% x = reshape(F, [],1);
% y = reshape(A, [],1);

%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri di controllo
%<<<<<<<<<<<<<<<<<<<<<<<<

% Parametri fisici
[g,div_F,div_A,fs] = parametri_fisici();

% Parametri di ricerca
[soglia,delay,inizio,fine] = parametri_ricerca_picchi(fs,x);

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
ascissamin=20;          % Frequenza minima da plottare nei grafici
ascissamax=20000;       % Frequenza massima da plottare nei grafici

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


[pks,picchi_t,w,p] = findpeaks(x,fs,'MinPeakProminence',20,'MinPeakDistance',0.3,'Annotate','extents');
picchi_t=picchi_t(w<0.2);
pks=pks(w<0.2);
picchi=picchi_t*fs;

figure(1)
subplot(2,1,1), hold on, plot(picchi_t,pks,'*')

n_picchi=length(picchi)
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

% Generazione finestra Accelerazione
M_win = L_pre; % M_win è l'ordine della finestratura, la lunghezza del fronte di salita + discesa
L_plateau = 0.0625; % L_plateau indica la lunghezza del plateau rispetto alla massima possibile
switch wintype
    case 'hann'
        curve = hann(M_win);
        plateau_1 = ones(round((r-M_win-1)*L_plateau), 1);
        plateau_0 = zeros((r-M_win-1) - round((r-M_win-1)*L_plateau), 1);
        win_A = [curve(1:end/2); plateau_1; curve(end/2:end); plateau_0];
    case 'rect' %finestratura rettangolare
        curve = ones(M_win);
        win_A = [curve(1:end/2);ones((r-2*M_win-1),1);curve(end/2:end)];
    case 'none'
        win_A = ones(r,1);
end
E_win=sum(win_A.^2)/length(win_A);
dt=1/fs; time1=1000*(0:dt:r/fs-dt);

L_win_F=L_pre+round(1*fs/1000)+M_win/2;
win_F=[ones(L_win_F-(M_win/2),1);curve(end/2:end-1);zeros(Lsample-L_win_F,1)];
win_F(1,1)=0;

% Plot delle finestre
figure(101),hold on 
subplot(2,2,2)
plot(time1,win_A*20*log10(max(max(A))));
subplot(2,2,1);
plot(time1,win_F*max(max(F)));
hold off

%<<<<<<<<<<<<<<<<<<<<<
% Analisi smorzamento
%<<<<<<<<<<<<<<<<<<<<<
A_reverse=cumsum(flip(A.^2));
figure (2), hold on, 
plot(10*log10(movmean(A.^2,20)))
figure(4), hold on, plot(10*log10(flip(A_reverse)))

% Finestratura
A_filt=A.*win_A;
F_filt=F.*win_F;

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
[PSD_F, f]= periodogram(F_filt, win_1, r, fs); %PSD Forza [N^2]
[PSD_A, f]= periodogram(A, win_A, r, fs); %PSD Accelerazione [g^2]
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
[r,c]=size(PSD_F);
tagli=[];
scarti=0;
for jj=1:(c-scarti)
    if (max(PSD_F(:,jj-scarti))< 0.5*max(max(PSD_F)))
        PSD_F(:,jj-scarti)=[];
        tagli=[tagli; jj];
        scarti=scarti+1;
    end
end
picchi_sel2 = picchi_sel2 - scarti
PSD_A(:,tagli)=[];

%<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo Dstiff totale
%<<<<<<<<<<<<<<<<<<<<<<<

% Calcolo tramite periodogram
PSD_Fav_pgram = mean(PSD_F, 2);
PSD_Aav_pgram = mean(PSD_A, 2);
PSD_Vav_pgram = PSD_Aav_pgram./(1i*2*pi*f).^2; % velocità
PSD_Dav_pgram = PSD_Vav_pgram./(1i*2*pi*f).^2; % displacement
PSD_Kav_pgram = PSD_Fav_pgram./PSD_Dav_pgram; % dynamic stiffness calcolata tramite periodogram

% Calcolo tramite FFT
FFT_Fav_fft = mean(FFT_F,2);
FFT_Aav_fft = mean(FFT_A,2);
FFT_Vav_fft = FFT_Aav_fft./(1i*2*pi*f_fft);
FFT_Dav_fft = FFT_Aav_fft./(1i*2*pi*f_fft).^2;
FFT_Kav_fft = FFT_Fav_fft./FFT_Dav_fft;

% Plot Forza media vs Accelerazione media e confronto tra periodogram e fft
i=1;
figure(110), hold on

subplot (2,1,1),hold on,
plot(f_fft,10*log10(mean(PSD_F_fft,2)),'color',string(colore(i,:)))
plot(f,10*log10(PSD_Fav_pgram(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})',...
    'PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

subplot (2,1,2),hold on
plot(f_fft,10*log10(mean(PSD_A_fft,2)),'color',string(colore(i,:)))
plot(f,10*log10(PSD_Aav_pgram(1:length(f),1)),'-.','color',string(colore(i,:)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
legend('PSD tramite Fft (2\cdotabs(Fft(F))^2/ length(Fft(F(:,1)))/fs/E_{win})',...
    'PSD tramite Periodogram')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

saveas(gcf,'Forza_vs_Accelerazione.fig');

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% COERENZA usando Forza / Accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%Coerenza usando tutto il segnale campionato e filtrato
F_filtall = reshape(F_filt, [],1);
A_filtall = reshape(A_filt, [],1);
[r,c]=size(F_filt);
[Cxy1,f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs);
save (cell2mat(['Coerenza_',conf.campione,'_',conf.piastra]), 'Cxy1');

PSD_Kav_gate=PSD_Kav_pgram;
f_Kav_gate=f;
tagli=[];
for ii=1:length(f)
    if Cxy1(ii)<0.7 & Cxy1(ii)>0.3
        tagli=[tagli,ii];
    end
end
PSD_Kav_gate(tagli)=[];
f_Kav_gate(tagli)=[];


Acc=PSD_Aav_pgram./PSD_Fav_pgram;
save(('Acc'), 'Acc');
save(('Forza'), 'PSD_Fav_pgram');
save(('Accelerazione'), 'PSD_Aav_pgram');

figure(111), hold on

subplot (2,1,1),hold on,
plot(f_coer,Cxy1,'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
legend('Coerenza')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

subplot (2,1,2),hold on
plot(f,Acc(1:length(f)),'LineWidth',2)
grid on, set(gca, 'XScale', 'log'), xlim([ascissamin ascissamax])
legend('Accelerazione su Forza')%,'2*fft^2/(l*fs)','2*fft^2/(l*fs)*E_win')

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo della f massima media
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
PSD_Fav_dB = 10*log10(PSD_Fav_pgram);
fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
fmax=f(fmax(1)+f0);
%plot sulla PSD della forza
figure (111)
subplot (2,1,2),hold on
%xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
hold off
saveas(gcf,'Acc.fig');

fr=f(find (Acc==max(Acc(find (f>ascissamin,1):find (f>fmax,1)))));

m=piastre.massa(conf.piastra); %massa della piastra in uso
h=campioni.h(conf.campione);
s=pi*(campioni.d(conf.campione)/2)^2;

K0_av=(2*pi*fr)^2*m
St=K0_av/s
E_av=K0_av*h/s

Risultati_Acc = table(fr,K0_av,St,E_av);
Risultati_Acc.Properties.VariableNames={'f_r','K_0','S_t','E'}
save (['Risultati_Acc.mat'], 'Risultati_Acc');

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Confronto PSD calcolata tramite Fft e Psd
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (109),hold on,
%sgtitle {'Dynamic Stiffness media K(f)'}

subplot(4,1,[2,3]),hold on
plot (f_fft,10*log10(abs(FFT_Kav_fft).^2),'color',string(colore(i,:)),'LineWidth',1)
plot (f(1:length(PSD_Kav_pgram)),10*log10(PSD_Kav_pgram),'.','color',string(colore(i,:)),'LineWidth',1)
semilogx (f_Kav_gate,10*log10(abs(PSD_Kav_gate).^2),'r.','LineWidth',2)
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

%<<<<<<<<<<<<<<<<<<<<<<<
% Plot del massimi di F
%>>>>>>>>>>>>>>>>>>>>>>>
PSD_F2=PSD_F;
PSD_A2=PSD_A;
% Ricerca dei massimi e ordinamento
PSD_Fsort=[];
PSD_Asort=[];
for i=1:picchi_sel2
    pos=find(sqrt(max(PSD_F2))==sqrt(max(max(PSD_F2))));
    PSD_Fsort=[PSD_Fsort,PSD_F2(:,pos)];
    PSD_Asort=[PSD_Asort,PSD_A2(:,pos)];
    PSD_F2(:,pos)=[];
    PSD_A2(:,pos)=[];
end
% plot dei massimi
figure (5), hold on
ylimsup =  floor(max(max(PSD_F)));
ylim=([0 1]); 
plot(max(PSD_Fsort));
plot(max(PSD_F));
hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi statistica post filtraggio
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[Y,E] = discretize(sqrt(max(abs(PSD_F))),bin);
values=1:bin;

%figure (3), hold on
%histfit(sqrt(max(abs(PSD_F))),bin);
%hold off

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
kkk=0;
[R,C]=size(PSD_F);

for indice = 1:bin
    kkk=kkk+1;

    PSD_Fbin=[];
    PSD_Abin=[];

    % Seleziono i campioni appartenenti al BIN
    for jj=1:C
        if Y(jj)==indice
            PSD_Fbin=[PSD_Fbin,PSD_F(:,jj)]; %#ok<AGROW>
            PSD_Abin=[PSD_Abin,PSD_A(:,jj)]; %#ok<AGROW>
        end
    end
    
    % Calcolo dimensioni del BIN
    [RR,CC]=size(PSD_Fbin);
    
    if CC>=1 % se il BIN non è vuoto
        
        F_bin=[F_bin,(E(kkk)+E(kkk+1))/2]; %#ok<AGROW> %Fbin è la forza media del BIN
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo la media degli spettri
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_bin = mean(PSD_Fbin,2);
        PSD_Aav_bin = mean(PSD_Abin,2);
        PSD_Vav_bin = PSD_Aav_bin./(1i*2*pi*f).^2; %velocità
        PSD_Dav_bin = PSD_Vav_bin./(1i*2*pi*f).^2; %spostamento
        % Dynamic Stiffness
        PSD_Kav_bin = PSD_Fav_bin./PSD_Dav_bin; %Modulo della Dynamic Stiffness
        save (['Dstiffness_av_bin.mat'], 'PSD_Kav_bin');
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 10*log10(PSD_Fav_bin);
        fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax=f(fmax(1)+f0);        
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot di segnali e spettri (PSD)
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure(101), grid on,
        %sgtitle(misura)
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
        % Plot spettri medi
        subplot(2,2,3), hold on,
        plot (f, 10*log10(PSD_Fav_bin), '-.','color',string(colore(kkk,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on,
        plot (f, 10*log10(PSD_Aav_bin), '-.' ,'color',string(colore(kkk,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        %plot frequenza massima sulla PSD della forza
        subplot (2,2,3), hold on
        %xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
        hold off
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plod Dstiff singole in grafico per BIN
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (round(indice+200))
        hold on
        % Ciclo per plottare le singole dstiff
        for iii=1:CC
            semilogx (f, 10*log10(abs(PSD_Fbin(:,iii).*((1i*2*pi*f).^4)./PSD_Abin (:,iii) )),...
                'color',string(colore(kkk,:)), 'LineWidth', 1),    
        end
        % Plot valore medio del bin
        semilogx (f,    10*log10(PSD_Kav_bin),'k--', 'LineWidth', 2),        
        % Plot frequenza massima
        %xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        % Settaggio parametri Grafico
        set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
        titolo=['Dynamic Stiffness (Force/Displacement) Amplitude (fron ',num2str(E(indice)),...
            ' to ',num2str(E(indice+1)),' N)'];
        title([titolo]),
        grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
        % Salvataggio
        saveas (gcf, ['from ',num2str(E(indice)),' to ',num2str(E(indice+1)),' N)',' N Dstiff.fig'])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot grafico riassuntivo con dstiff medi di tutti i BIN
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (107),hold on,
        semilogx (f, 10*log10(PSD_Kav_bin),'color',string(colore(kkk,:)), 'LineWidth', 1),
        %semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3),
        %xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        
        %<<<<<<<<<<<<<<<<<<<<
        % Calcolo K0 del bin
        %<<<<<<<<<<<<<<<<<<<<
        min_k = min(PSD_Kav_bin(lim_inf(end):lim_sup(1)));
        fr = f(lim_inf(end)+find(PSD_Kav_bin(lim_inf(end):lim_sup(1)) == min_k));
        K0_av_bin=[K0_av_bin,(2*pi*fr)^2*m];
        E_av_bin=[E_av_bin,K0_av_bin*h/s];
        PSD_K_bin = [PSD_K_bin,K0_av_bin];

    end
    
end

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
figure (107),  hold on
plot (f_fft, 10*log10(abs(FFT_Kav_fft).^2),'LineWidth',2)
plot (f(1:length(PSD_Kav_pgram)), 10*log10(PSD_Kav_pgram),'.','LineWidth',2)
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
xlabel('Frequenza [Hz]'), ylabel('Psd Rigidezza Dinamica [dB ref 1 N/m]'),
title(['Pds della rigidezza dinamica [N/mHz]. Campione: ',cell2mat([conf.campione,...
    ' + ',conf.piastra,'. Adesivo: ',conf.adesivo,'.'])]),
% Salvataggio
saveas (gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.fig']))
saveas (gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.png']))

% Salvataggio delle dstiffness medie FFT e PSD
save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Psd.mat']),'PSD_Kav_pgram');
save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Fft.mat']),'FFT_Kav_fft');


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

%<<<<<<<<<<<<<<<<<<<<<<<<<
% Stampo K(F) interpolato
%<<<<<<<<<<<<<<<<<<<<<<<<<
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

%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% test della crossrasformata per la media delle risposte
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

[pxy, ~] = cpsd (sng(1).F_filt, sng(1).A_filt, [], [], f, fs);

% pxx = cpsd (sng(1).F_filt, sng(1).F_filt, [], [], f, fs);
pxx = PSD(1).F;

figure (300)
hold on
plot (f, 20*log10(pxy(:,1)./pxx(:,1)))
plot (f, 10*log10(pxx(:,1)./PSD.))
set(gca, 'XScale', 'log')