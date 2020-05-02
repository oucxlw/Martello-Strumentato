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
appoggio={'cemento'};
adesivo={'gesso'};
punta={'metallo'};
martellatore=2;
accelerometro=0;
conf = table(campione, piastra, appoggio, adesivo, punta, martellatore, accelerometro)

%%
clc
close all
clear variables

dati = load ('dati1');
% dati(2) = load ('dati2');
% dati(3) = load ('dati3');
% dati(4) = load ('dati1a');
% dati(5) = load ('dati2a');
% dati(6) = load ('dati3a');
% dati(7) = load ('dati1b');
% dati(8) = load ('dati2b');
% dati(9) = load ('dati3b');
% dati(10) = load ('dati1c');
% dati(11) = load ('dati2c');
% dati(12) = load ('dati3c');
conf = dati(1).conf;
fs=dati(1).fs;

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf,piastre);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
misure = cell(length(dati),1);
sng = struct('x',misure,'y',[]);
N = length (dati);
for i = 1:N
    sng(i).x = dati(i).data(:,1);
    sng(i).y = dati(i).data(:,conf.accelerometro+2);
end

%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri di controllo
%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri fisici
[g,div_F,div_A,fs] = parametri_fisici();

% Parametri di ricerca
fine = cell(1,N);
for i = 1:length(sng)
    [soglia,delay,inizio,fine{i,1}] = parametri_ricerca_picchi(fs,sng(i).x);
end

% Parametri di creazione dei campioni
[Lsample,L_pre,L_coda] = parametri_creamatrice();

% Parametri di filtro
bandwidth = 0;
filt_doppi = 0;           % Se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi

% Finestratura
wintype = 'hann';
% Determina il tipo di finestratura da utilizzare per i segnali di
% accelerazione:
% hann = utilizza la finestra di hann;
% rect = utilizza una finestratura rettangolare;
% none = non applica nessuna finestratura.

% Plotting
ascissamin = 20;          % Frequenza minima da plottare nei grafici
ascissamax = 5000;       % Frequenza massima da plottare nei grafici

misura = cell2mat(['Campione ',conf.campione,', martellatore ',...
    num2str(conf.martellatore),', punta ',conf.punta,', piastra ',conf.piastra,...
    ',Hann,PSDvsFFT, ',num2str(bandwidth),' Hz']);

colore = [
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
for i = 1:N
    sng(i).x = div_F * sng(i).x;
    sng(i).y = g * div_A * sng(i).y;
end

% Plot dei segnali
for i=1:N
    L = length(sng(i).x);
    dt=1/fs; time=[0:dt:L/fs-dt];
    figure(i)
    subplot(2,1,1), hold on, plot (time, sng(i).x)
    subplot(2,1,2), hold on, plot (time, sng(i).y),
    hold off
end

%<<<<<<<<<<<<<<<<<<<<
% Ricerca dei PICCHI
%<<<<<<<<<<<<<<<<<<<<
prominanza=10;%25;
distanza=0.5;
larghezza=10;
fpass=35
x_Hpass = cell(1,N);
for i=1:N
    x_Hpass{i} = highpass(sng(i).x,fpass,fs);
end

picchi_t = cell(1,N);
picchi1 = cell(1,N);
p = cell(1,N);
w = cell(1,N);
n_picchi = cell(1,N);

for i=1:N
    [pks1,picchi_t{i},w{i},p{i}] = findpeaks(x_Hpass{i},fs,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
    
    picchi_t{i}=picchi_t{i}(w{i}<larghezza);
    pks1=pks1(w{i}<larghezza);
    picchi1{i} = picchi_t{i}*fs;
    n_picchi{i} = length(picchi1{i})
    
    figure(i)
    subplot(2,1,1), hold on, plot(picchi_t{i},pks1,'*')
    findpeaks(x_Hpass{i},fs,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
end
%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pos = cell(1:N);
picchi_sel1 = zeros(1,N);
for i = 1:N
    [sng(i).F, pos{i}] = creamatriceforza_noavg (sng(i).x, picchi1{i}, n_picchi{i}, L_pre, L_coda, filt_doppi, fs);
    [sng(i).A] = creamatriceaccelerazione (sng(i).y, pos{i}, L_pre, L_coda, fs);
    picchi_sel1(i) = length(pos{i});
end
picchi_sel1



% % Accodo tutti i segnali di forza e accelerazione in un'unica matrice
% F=[F1 F2 F3];
% A=[A1 A2 A3];

% Matrice degli indici
% Se l'elemento i-esimo di I vale x, la colonna i-esima di F appartiene
% alla matrice forza Fx.
% I=[ones(size(F1(1,:))) 2*ones(size(F2(1,:))) 3*ones(size(F3(1,:)))];
% I0=I;

[r,c] = size(sng(1).F(:,1));
%<<<<<<<<<<<<<<
% Finestratura
%<<<<<<<<<<<<<<
% Generazione finestra Accelerazione
M_win = L_pre; %M_win � l'ordine della finestratura, ossia la sua lunghezza
switch wintype
    case 'hann'
        curve=hann(M_win);
        win_A = [curve(1:end/2);ones((r-M_win-1),1);curve(end/2:end)];
    case 'rect' %finestratura rettangolare
        curve = ones(M_win);
        win_A = [curve(1:end/2);ones((r-2*M_win-1),1);curve(end/2:end)];
    case 'none'
        win_A = ones(r,1);
end
E_win=sum(win_A.^2)/length(win_A);
dt = 1/fs; time1 = 1000*(0:dt:r/fs-dt);

% Generazione finestra Forza
L_win_F = L_pre+round(1*fs/1000)+M_win/2;
win_F = [ones(L_win_F-(M_win/2),1);curve(end/2:end-1);zeros(Lsample-L_win_F,1)];
win_F(1,1) = 0;

% Plot delle finestre
figure(101),hold on
subplot(2,2,2)
plot(time1, win_A * 20*log10(max(max(sng(1).A))));
subplot(2,2,1);
plot(time1, win_F * max(max(sng(1).F)));
hold off

% Finestratura
for i = 1:N
    sng(i).F_filt = sng(i).F  .* win_F;
    sng(i).A_filt = sng(i).A  .* win_A;
end

%<<<<<<<<<<<<<<<<<<<<<
% Analisi smorzamento
%<<<<<<<<<<<<<<<<<<<<<
% A_reverse=cumsum(flip(A.^2));
% figure (14), hold on,
% plot(10*log10(movmean(A.^2,20)))
% figure(15), hold on, plot(10*log10(flip(A_reverse)))

%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle PSDs
%<<<<<<<<<<<<<<<<<<<<
% Trasformata (non normalizzata rispetto al numero di punti fft)
FFT = struct('F',cell(1,N),'A',[]);
for i = 1:N
    FFT(i).F = fft(sng(i).F);
    FFT(i).A = fft(sng(i).A);
end

% Frequenza della FFT_F
f_fft = 0:(r-1);
f_fft = f_fft'/(r-1)*fs;


% Calcolo della PSD tramite FFT
PSD = struct('F',cell(1,N),'A',[],'F_fft',[],'A_fft',[]);

for i = 1:N
    PSD(i).F_fft = abs(FFT(i).F).^2./length (FFT(i).F(:,1))/fs/1; %PSD di F stimata tramite fft
    PSD(i).F_fft(2:end,:) = 2*PSD(i).F_fft(2:end,:);
    PSD(i).A_fft = abs(FFT(i).A).^2./length (FFT(i).A(:,1))/fs/E_win;
    PSD(i).A_fft(2:end,:) = 2*PSD(i).A_fft(2:end,:);
end

% Calcolo della PSD tramite Periodogram
win_1=ones(size(win_F)); % finestra unitaria per non far calcolare la normalizzazione a periodogram

% Periodogram non gestisce pi� di 40 colonne contemporaneamente quindi
% eseguo le operazioni ciclicamente
for i=1:N
    [PSD(i).F, ~]= periodogram(sng(i).F_filt, win_1, r, fs); %PSD Forza [N^2]
    [PSD(i).A, f]= periodogram(sng(i).A, win_A, r, fs); %PSD Accelerazione [g^2]
end

%<<<<<<<<<<<<<<<<<<<<<<
% Filtraggio per Banda
%<<<<<<<<<<<<<<<<<<<<<<
f0=find(f>50,1); % posizione che corrisponde a 50 Hz
scarti = zeros(1,N);
if bandwidth>0
    tagli = [];
    for i=1:N % Scorro sulle misure
        [~,C]=size(PSD(i).F);
        for jj = 1:C % Scorro sui colpi
            fmax = find(PSD(i).F(f0:end, jj) < ((PSD(i).F(f0, jj)/10)), 1);
            fmax = f(fmax+f0);
            if  fmax < bandwidth
                tagli = [tagli, jj]; %#ok<AGROW>
            end
            PSD(i).F(:, tagli) = [];
            PSD(i).A(:, tagli) = [];
            PSD(i).F_fft(:, tagli) = [];
            PSD(i).A_fft(:, tagli) = [];
            sng(i).F(:, tagli) = [];
            sng(i).A(:, tagli) = [];
            
            sng(i).F_filt(:, tagli) = [];
            sng(i).A_filt(:, tagli) = [];
            
            scarti(i) = length(tagli);
        end
    end
end
picchi_sel2 = picchi_sel1 - scarti

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
% Filtro i picchi che sono molto pi� piccoli del picco pi� grande
scarti = zeros(1,N);
for i = 1:N
    tagli=[];
    [~,C] = size(PSD(i).F);
    for jj = 1:C
        if (PSD(i).F(f0,jj) < 0 * max(PSD(i).F(f0,:)))
            tagli = [tagli, jj]; %#ok<AGROW>
        end
    end
    PSD(i).F(:, tagli) = [];
    PSD(i).A(:, tagli) = [];
    PSD(i).F_fft(:, tagli) = [];
    PSD(i).A_fft(:, tagli) = [];
    sng(i).F(:, tagli) = [];
    sng(i).A(:, tagli) = [];
    sng(i).F_filt(:, tagli) = [];
    sng(i).A_filt(:, tagli) = [];
    scarti(i) = length(tagli);
    
end
picchi_sel2 = picchi_sel1 - scarti

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo delle grandezze fisiche delle singole registrazioni
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for i = 1:N
    % psd
    PSD(i).Fav = mean(PSD(i).F, 2);
    PSD(i).Aav = mean(PSD(i).A, 2);
    
    PSD(i).Vav = PSD(i).Aav./(1i*2*pi*f).^2; % velocita'
    PSD(i).Dav = PSD(i).Vav./(1i*2*pi*f).^2; % displacement
    
    PSD(i).Kav = PSD(i).Fav./PSD(i).Dav; % dynamic stiffness calcolata tramite periodogram
    
    % fft
    FFT(i).Fav = mean(FFT(i).F, 2);
    FFT(i).Aav = mean(FFT(i).A, 2);
    
    FFT(i).Vav = FFT(i).Aav ()./(1i*2*pi*f_fft); % velocita'
    FFT(i).Dav = FFT(i).Aav./(1i*2*pi*f_fft).^2; % displacement
    
    FFT(i).Kav = FFT(i).Fav ./ FFT(i).Dav;  % dynamic stiffness calcolata tramite fft
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo risultati - Dstiff totale
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo tramite periodogram
pnt_misura.PSD = struct('F',[],'A',[]);

pnt_misura.PSD.F = mean(vertcat([PSD([1:N]).F]),2);
pnt_misura.PSD.A = mean(vertcat([PSD([1:N]).A]),2);

pnt_misura.PSD.V = pnt_misura.PSD.A./(1i*2*pi*f).^2; % velocità
pnt_misura.PSD.D = pnt_misura.PSD.V./(1i*2*pi*f).^2; % displacement

pnt_misura.PSD.K = pnt_misura.PSD.F./pnt_misura.PSD.D; % Dynamic stiffness FFT
pnt_misura.PSD.MI = pnt_misura.PSD.F./pnt_misura.PSD.V; % Mechanical Impedence PSD

% Calcolo tramite fft
pnt_misura.FFT = struct('F',[],'A',[]);

pnt_misura.FFT.F = mean(vertcat([FFT([1:N]).F]),2);
pnt_misura.FFT.A = mean(vertcat([FFT([1:N]).A]),2);

pnt_misura.FFT.V = pnt_misura.FFT.A./(1i*2*pi*f_fft).^2; % velocità
pnt_misura.FFT.D = pnt_misura.FFT.V./(1i*2*pi*f_fft).^2; % displacement

pnt_misura.FFT.K = pnt_misura.FFT.F./pnt_misura.FFT.D; % Dynamic stiffness FFT
pnt_misura.FFT.MI = pnt_misura.FFT.F./pnt_misura.FFT.V; % Mechanical Impedence PSD


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

%<<<<<<<<<<<<<<<<<<<<<
% Analisi bin per bin
%<<<<<<<<<<<<<<<<<<<<<
m = piastre.massa(conf.piastra); %massa della piastra in uso
h = campioni.h(conf.campione);
s = pi*(campioni.d(conf.campione)/2)^2;


% [R,C]=size(PSD_F);
F_max = [];
fmax_misura = []; % la collzione della frequenza massima della bandwidth della forza a -10dB

result.indici_colpo.max_A = struct('fd', misure, 'K0', misure,'S_star', [], 'E', []);
result.indici_colpo.min_K = struct('fd', misure, 'K0', misure,'S_star', [], 'E', []);
result.indici_colpo.UNI_10570_1997 = struct('fd', misure, 'nu_d', [], 'K0', [],'S_star', [], 'E', []);

% Struttura con i risultati sintetici di ogni singola misura 
% Effettuati sui segnali mediati per sciascuna singola registrazione
result.indici.max_A = struct('fd', misure, 'K0', misure,'S_star', [], 'E', []);
result.indici.min_K = struct('fd', misure, 'K0', misure,'S_star', [], 'E', []);
result.indici.UNI_10570_1997 = struct('fd', misure, 'nu_d', [], 'K0', [],'S_star', [], 'E', []);

% Struttura con i risultati sintetici di tutto l'insieme di misure
% effettuato sul segnale ottenuto mediando tutte le registrazioni
pnt_misura.indici.max_A = struct('fd', misure, 'K0', [],'S_star', [], 'E', []);
pnt_misura.indici.min_K = struct('fd', misure, 'K0', [],'S_star', [], 'E', []);
pnt_misura.indici.UNI_10570_1997 = struct('fd', [], 'nu_d', [], 'K0', [],'S_star', [], 'E', []);

% Parametri di ingresso:
% PSD_F, PSD_A, F_filt, A_filt,I
for mis = 1:N
    % Calcolo dimensioni del BIN
    [RR,CC]=size(PSD(mis).F);
    
    if CC >= 1 % se il BIN non � vuoto
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo Forza media di ciascuna registrazione
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sng(mis).F_max_av = mean(max(sng(mis).F));
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % COERENZA usando Forza / Accelerazione
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sng(mis).F_filtall = reshape(sng(mis).F_filt, [],1);
        sng(mis).A_filtall = reshape(sng(mis).A_filt, [],1);
        [~,c]=size(sng(mis).F_filt);
        
        % Calcolo coerenza
        [sng(mis).Cxy,~] = mscohere(sng(mis).F_filtall, sng(mis).A_filtall, round(length(sng(mis).F_filtall)./c),[],f_fft,fs);
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo deviazione standard
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD(mis).devst_F = std (abs(PSD(mis).F), 0, 2);
        PSD(mis).devst_A = std (abs(PSD(mis).A), 0, 2);
        
        dF = PSD(mis).devst_F;
        dA = PSD(mis).devst_A;
        b  = abs((1i*2*pi*f).^4);
        FF = abs(PSD(mis).Fav);
        AA = abs(PSD(mis).Fav);
        C  = abs(sng(mis).Cxy(1:1+end/2));
        KK = PSD(mis).Kav;
        
        dK = ( ((b./AA).*dF).^2 + ((-b.*FF./AA.^2).*dA).^2 + C.*(b./AA).*(-b.*FF./AA.^2).*dF.*dA  ).^(1/2);
        PSD(mis).devst_K = dK;
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima delle forza media
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 10*log10(PSD(mis).Fav);
        fmax = find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax = f(fmax(1)+f0);
        PSD(mis).fmax = fmax;
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot di segnali e spettri (PSD)
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure(101), grid on,
        sgtitle(misura)
        for i = 1:c
            subplot(2,2,1), hold on
            plot (time1, sng(mis).F(:,i),'color',string(colore(mis,:))),
            
            subplot(2,2,2), hold on,
            plot (time1, 20*log10(abs(sng(mis).A(:, i))),'color',string(colore(mis,:)))
            
            subplot(2,2,3), hold on
            semilogx (f, 10*log10(abs(PSD(mis).F(:, i))),'--','color',string(colore(mis,:)))
            
            subplot(2,2,4), hold on,
            semilogx (f, 10*log10(abs(PSD(mis).A(:, i))),'--','color',string(colore(mis,:)))
        end
        
        % Plot spettri medispo
        subplot(2,2,3), hold on,
        plot (f, 10*log10(abs(PSD(mis).Fav)), 'color', string(colore(mis,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on,
        plot (f, 10*log10(abs(PSD(mis).Aav)), 'color', string(colore(mis,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        %plot frequenza massima sulla PSD della forza
        subplot (2,2,3), hold on
        %xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
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
        semilogx (f, 10 * log10(abs(PSD(mis).Kav)),'color',string(colore(mis,:)), 'LineWidth', 2),
        % Plot deviazione standard
        semilogx (f, 10 * log10(abs(PSD(mis).Kav) + dK), '-.','color',string(colore(mis,:)), 'LineWidth', 1),
        semilogx (f, 10 * real(log10(abs(PSD(mis).Kav) - dK)), '-.','color',string(colore(mis,:)), 'LineWidth', 1),
        
        % Plot limite in frequenza
        xline(PSD(mis).fmax, '.', ['Limite in frequenza: ',num2str(round(PSD(mis).fmax)),' Hz'],'color',string(colore(mis,:)));
        
        % Plot fase
        subplot(4,1,4), hold on
        semilogx (f_fft, 180/pi * unwrap(angle(FFT(mis).Kav)), 'color',string(colore(mis,:)),'LineWidth', 2)
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot Accelerazione e forza
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure(110), hold on
        subplot (2,1,1),hold on,
        plot (f, 10 * log10(FF),'color',string(colore(mis,:)),'LineWidth',2)
        plot (f, 10 * log10(FF + dF),'--','color',string(colore(mis,:)),'LineWidth',1)
        plot (f, 10 * real(log10(FF - dF)),'--','color',string(colore(mis,:)),'LineWidth',1)
        
        subplot (2,1,2),hold on,
        plot (f, 10 * log10(AA),'color',string(colore(mis,:)),'LineWidth',2)
        plot (f, 10 * log10(AA + dA),'--','color',string(colore(mis,:)),'LineWidth',1)
        plot (f, 10 * real(log10(AA - dA)),'--','color',string(colore(mis,:)),'LineWidth',1)
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo K0 della misura
        %<<<<<<<<<<<<<<<<<<<<<<<<<
        % sono arrivato qui
        
        lim_sup =  find(f > min(vertcat([PSD.fmax]))*0.35);
        lim_inf =  find(f < 20);
        
        % secondo minimo di K
        min_k = min(PSD(mis).Kav(lim_inf(end):lim_sup(1)));
        fd = f( lim_inf(end)+find(PSD(mis).Kav(lim_inf(end):lim_sup(1)) == min_k) );
        result.indici.min_K(mis).fd = fd;
        result.indici.min_K(mis).K0 = (2*pi*fd)^2 * m;
        result.indici.min_K(mis).S_star = (2*pi*fd)^2 * m/s;
        result.indici.min_K(mis).E  = (2*pi*fd)^2 * m * h/s;
        
        % seconfo massimo di A
        temp = PSD(mis).Aav./PSD(mis).Fav;
        max_k = max(abs(temp(lim_inf(end):lim_sup(1))));
        fd = f(lim_inf(end)+find(temp(lim_inf(end):lim_sup(1)) == max_k) );
        result.indici.max_A(mis).fd = fd;
        result.indici.max_A(mis).K0 = (2*pi*fd)^2 * m;
        result.indici.max_A(mis).S_star = (2*pi*fd)^2 * m/s;
        result.indici.max_A(mis).E  = (2*pi*fd)^2 * m * h/s;
    end
end

%%

figure(110)
saveas(gcf,'Forza_vs_Accelerazione.fig');

% Salvataggio coerenza
% save (cell2mat(['Collezione_Coerenza_',conf.campione,'_',conf.piastra]), 'Cxy');

% % Salvataggio rigidezza dinamica PSD
% save (cell2mat(['Collezione_Dstiffness_PSD_',conf.campione,'_',conf.piastra]), 'PSD_Kav_misura');
% 
% % Salvataggio rigidezza dinamica FFT
% save (cell2mat(['Collezione_Dstiffness_FFT_',conf.campione,'_',conf.piastra]), 'FFT_K_misura');
% 
% % Salvataggio deviazione standard rigidezza dinamica PSD
% save (cell2mat(['Collezione_DevSt_',conf.campione,'_',conf.piastra]), 'devst_K_misura');


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo delle grandezze sintetiche della misura
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% 
% F_filtall = reshape(F_filt, [],1);
% A_filtall = reshape(A_filt, [],1);
% [r,c]=size(F_filt);
% 
% % Calcolo coerenza
% [Cxy_sintetico,f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs);

% deviazione standard
pnt_misura.PSD.devst_K = mean(vertcat([PSD([1:N]).devst_K]).^2, 2);
pnt_misura.PSD.devst_F = std (vertcat([PSD([1:N]).F]),0,2);
pnt_misura.PSD.devst_A = std (vertcat([PSD([1:N]).A]),0,2);

% dF = pnt_misura.PSD.devst_F;
% dA  = pnt_misura.PSD.devst_A;
% b = abs((1i*2*pi*f).^4);
% FF = abs(mean(PSD_F,2));
% AA = abs(mean(PSD_A,2));
% C = abs(Cxy_sintetico(1:1+end/2));
% KK = PSD_Kav_misura(:,mis);
% % K = F*(1i*2*pi).^4 /A
% 
% dK = ( ((b./AA).*dF).^2 + ((-b.*FF./AA.^2).*dA).^2 + C.*(b./AA).*(-b.*FF./AA.^2).*dF.*dA  ).^(1/2);
% devst_K_sintetico = dK;
% % Salvataggio deviazione standard sintetica rigidezza dinamica PSD
% save (cell2mat(['DevSt_sintetica_',conf.campione,'_',conf.piastra]), 'devst_K_sintetico');

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

% % Salvataggio delle k0 misurate bin per bin
% save (cell2mat(['Dstiffness_bin_',conf.campione,'_',conf.piastra,'_Fft.mat']),'PSD_K_bin');

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dstiff totale e settaggio parametri grafico
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (107)
% sgtitle(['PSD della rigidezza dinamica [N/mHz]. Campione: ',cell2mat([conf.campione,...
%     ' + ',conf.piastra,'. Adesivo: ',conf.adesivo,'.'])])

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
% legend(['misura 1 (',num2str(round(F_max_av (1))),' N)'],...
%     ['misura 2 (',num2str(round(F_max_av (2))),' N)'],...
%     ['misura 3 (',num2str(round(F_max_av (3))),' N)'])

% Salvataggio
saveas(gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.fig']))
%saveas(gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.png']))

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Accelerazione e forza
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure(110), hold on

subplot(2,1,1),hold on
set(gca, 'XScale', 'log')%, set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Forza [dB @ 1 N]'),

subplot(2,1,2),hold on,
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Accelerazione [dB @ 1 m/s^2]'),
saveas(gcf, cell2mat(['Forza_vs_Accelerazione_',conf.campione,'_',conf.piastra,'.fig']))


% Salvataggio delle dstiffness medie FFT e PSD
% save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Psd.mat']),'PSD_Kav_pgram');
% save (cell2mat(['Dstiffness_',conf.campione,'_',conf.piastra,'_Fft.mat']),'FFT_Kav_fft');

figure (106)
grid on %xlim([0 10])
set(gca, 'XScale', 'log')
xlim([10 20000])
hold on

for i=1:3
    subplot(2,1,1),hold on
    title('Impedenza meccanica')
    plot(f,10*log10(abs(PSD_Fav_misura(:,i)./PSD_Vav_misura(:,i))))
    subplot(2,1,2),hold on
    plot(f_fft,180/pi*unwrap(angle(FFT_F_misura(:,i)./FFT_V_misura(:,i))))
    
end
subplot(2,1,1),hold on, set(gca, 'XScale', 'log')

subplot(2,1,2),hold on, set(gca, 'XScale', 'log')
ylim([-200 200])


% accelerazione/forza
f1=100
f2=2000
f_zoom = f1:0.5:f2;
PSD_F_zoom = [];
PSD_A_zoom = [];

for i=1:3
    [PSD_Ftemp, f2,pxxc] = periodogram(F_filt(:,I==i), win_1, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Forza [N^2]
    [PSD_Atemp, f2,pyyc] = periodogram(A(:,I==i), win_A, f_zoom, fs, 'ConfidenceLevel',0.60); %PSD Accelerazione [g^2]
    PSD_F_zoom(:,i) = mean(PSD_Ftemp, 2);
    PSD_A_zoom(:,i) = mean(PSD_Atemp, 2);
end


figure (107)
subplot(4,1,[2,3]), hold on,
for i=1:3
    plot (f_zoom, 10*log10(abs(PSD_F_zoom(:,i)./PSD_A_zoom(:,i).*(1i*2*pi*f_zoom').^4)),'LineWidth', 3)
end
sgtitle(['PSD della rigidezza dinamica [N/mHz]. Campione: ',cell2mat([conf.campione,...
    ' + ',conf.piastra,'. Adesivo: ',conf.adesivo,'.'])])

hold off

saveas(gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.fig']))
saveas(gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.png']))

figure, hold on
plot(f,10*log10(PSD_Aav_misura(:,1)./PSD_Fav_misura(:,1)))
plot(f,10*log10(PSD_Aav_misura(:,2)./PSD_Fav_misura(:,2)))
plot(f,10*log10(PSD_Aav_misura(:,3)./PSD_Fav_misura(:,3)))
% psd zoom
plot(f_zoom, 10*log10(PSD_A_zoom(:,1)./PSD_F_zoom(:,1)))
plot(f_zoom, 10*log10(PSD_A_zoom(:,2)./PSD_F_zoom(:,2)))
plot(f_zoom, 10*log10(PSD_A_zoom(:,3)./PSD_F_zoom(:,3)))
plot(f_zoom, 10*log10(pyyc./PSD_F_zoom(:,3)), '-*')

grid on
set(gca, 'XScale', 'log'),
xlim([10 2000])
title('Accelerazione normalizzata','FontSize',18)
xlabel('Frequency [Hz]','FontSize',12),
ylabel('Acceleration/Force [Kg^{-1}]','FontSize',12),
saveas(gcf,'Accelerazione.fig');

% accelerazione/forza
figure (5), hold on
for i=1:3
    plot(f,PSD_Aav_misura(:,i)./PSD_Fav_misura(:,i))
end
grid on
set(gca, 'XScale', 'log')
xlim([10 1000])

% % velocit�/forza
% figure, hold on
% plot(f,PSD_Vav_misura(:,1)./PSD_Fav_misura(:,1))
% plot(f,PSD_Vav_misura(:,2)./PSD_Fav_misura(:,2))
% plot(f,PSD_Vav_misura(:,3)./PSD_Fav_misura(:,3))
% grid on
% set(gca, 'XScale', 'log')
% xlim([10 2000])
%
% % spostamento/forza
% figure, hold on
% plot(f,PSD_Dav_misura(:,1)./PSD_Fav_misura(:,1))
% plot(f,PSD_Dav_misura(:,2)./PSD_Fav_misura(:,2))
% plot(f,PSD_Dav_misura(:,3)./PSD_Fav_misura(:,3))
% grid on
% set(gca, 'XScale', 'log')
% xlim([10 2000])


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo frequenza di risonanza e K
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
m=piastre.massa(conf.piastra); %massa della piastra in uso
h=campioni.h(conf.campione);
s=pi*(campioni.d(conf.campione)/2)^2;
lim_sup =  find(f > min(vertcat([PSD.fmax]))*0.7);
lim_inf =  find(f < 100);
min_k = min(PSD_Kav_pgram(lim_inf(end):lim_sup(1)));

% basato su minimo Dynamic Stiffness
fd = f( lim_inf(end)+find(PSD_Kav_pgram(lim_inf(end):lim_sup(1)) == min_k)); % calcolata con la fft meno
K0_av = (2*pi*fd)^2*m
S_star_av = K0_av/s
E_av = K0_av*h/s

%basato sul massimo dell'accelerazione/forza
for i=1:3
    temp = PSD_Aav_misura(:,i)./PSD_Fav_misura(:,i);
    max_A = max(temp(lim_inf(end):lim_sup(1)));
    fr_A(i) = f(lim_inf(end)+find(temp(lim_inf(end):lim_sup(1)) == max_A)); % calcolata con la fft meno
end
K0_av_A = (2*pi*mean(fr_A))^2*m
S_star_av_A = K0_av_A/s
E_av_A = K0_av_A*h/s

figure (5), hold on
xl=xline (f(lim_sup(1)),'.',['limite superiore']); xl.LabelVerticalAlignment = 'bottom';
xl=xline (f(lim_inf(end)),'.',['limite superiore']); xl.LabelVerticalAlignment = 'bottom';

%figure(107),hold on, plot(f,20*log10(abs(K(1:length(f)))),'r.','LineWidth', 2)

%<<<<<<<<<<<<<<<<<<<<<<<<<
% Salvataggio dati misura
%<<<<<<<<<<<<<<<<<<<<<<<<<

save (['misura'],'Cxy','PSD_Kav_misura','FFT_K_misura','devst_K_misura','devst_K_sintetico',...
    'S_av_bin','conf','piastre','campioni','S_star_av_A')

save('frequenze.mat','f','f_fft')

save (cell2mat(['Dstiffness_bin_',conf.campione,'_',conf.piastra,'_Fft.mat']),'PSD_K_bin');


%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi risultati singoli colpi
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% accelerazione/forza
f1 = f(lim_inf(end));
f2 = f(lim_sup(1));
f_zoom = f1:0.2:f2;
prominanza=1;
distanza=10;

for i = 1:N
    [~,cc] = size(PSD(i).F);
    [PSD(i).F_zoom, ~, pxxc] = periodogram(sng(i).F_filt, win_1, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Forza [N^2]
    [PSD(i).A_zoom, ~, pyyc] = periodogram(sng(i).A, win_A, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Accelerazione [g^2]
    
    figure(120+i), hold on
    for j = 1:cc
        % secondo minimo di K
        PSD(i).K(:,j) = PSD(i).F(:,j) ./ PSD(i).A(:,j) .* (1i*2*pi*f).^4;
        plot (f, 10*log10(PSD(i).K(:,j)), 'color', string(colore(i,:)));
        
        PSD(i).K_zoom(:,j) = PSD(i).F_zoom(:,j) ./ PSD(i).A_zoom(:,j) .* (1i*2*pi*f_zoom').^4;
        plot (f_zoom, 10*log10(PSD(i).K_zoom(:,j)), 'color', string(colore(i+1,:)));
        
        min_k = min(PSD(i).K_zoom(lim_inf(end):lim_sup(1), j));
        fd = f_zoom(lim_inf(end) + find(PSD(i).K_zoom(lim_inf(end):lim_sup(1), j) == min_k) );
        result.indici_colpo.min_K(i).fd(j) = fd;
        result.indici_colpo.min_K(i).K0(j) = (2*pi*fd)^2 * m;
        result.indici_colpo.min_K(i).S_star(j) = (2*pi*fd)^2 * m/s;
        result.indici_colpo.min_K(i).E(j)  = (2*pi*fd)^2 * m * h/s;
    end
    xl=xline(f(lim_inf(end)),'.',['Lim inf ',f(lim_inf(end)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    xl=xline(f(lim_sup(1)),'.',['Lim sup ',f(lim_sup(1)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    grid on
    set(gca, 'XScale', 'log')
    xlim([ascissamin ascissamax])

    % secondo massimo di A
    figure(130+i), hold on
    for j = 1:cc
        % seconfo massimo di A
        temp = PSD(i).A(:,j) ./ PSD(i).F(:,j);
        plot (f, temp, 'color', string(colore(i,:)));
        
        temp_zoom = PSD(i).A_zoom(:,j) ./ PSD(i).F_zoom(:,j);
        plot (f_zoom, temp_zoom, 'color', string(colore(i+1,:)));
%                 
%         [pks,locs,wi,pi] = findpeaks(temp_zoom, f_zoom,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
%         plot(locs, pks,'*')

        max_A = max(temp_zoom);
        fd = f_zoom(temp_zoom == max_A );
        result.indici_colpo.max_A(i).fd(j) = fd;
        result.indici_colpo.max_A(i).K0(j) = (2*pi*fd)^2 * m;
        result.indici_colpo.max_A(i).S_star(j) = (2*pi*fd)^2 * m/s;
        result.indici_colpo.max_A(i).E(j)  = (2*pi*fd)^2 * m * h/s;
    end
    xl=xline(f(lim_inf(end)),'.',['Lim inf ',f(lim_inf(end)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    xl=xline(f(lim_sup(1)),'.',['Lim sup ',f(lim_sup(1)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    grid on
    set(gca, 'XScale', 'log')
    xlim([ascissamin ascissamax]), ylim([0 2*max_A])

end

%%
figure (120), hold on,
for i = 1:N
    [~,cc] = size(PSD(i).F);
    for j = 1:cc
        %plot (max(sng(i).F_filt(:,j)), result.indici_colpo.min_K(i).S_star(j)/10^6, '*', 'color', string(colore(i,:)))
        plot (max(sng(i).F_filt(:,j)), result.indici_colpo.max_A(i).S_star(j)/10^6, 'o', 'color', string(colore(i+1,:)), 'lineWidth', 3)
    end
end
grid on, %ylim([140 200])
X = max(vertcat([ sng([1:N]).F_filt]));
Y = vertcat([result.indici_colpo.max_A([1:N]).S_star]) / 10^6;

modelfun = @(b,x)b(1)*exp(-x/b(2))+b(3)
mdl = fitnlm(X, Y, modelfun, [100 50 10]);
[fit,R] = nlinfit(X(X>100 & X<400), Y(X>100 & X<400), modelfun, [200 150 150]);

X1 = 1:max(X*1.5);
plot (X1, feval(mdl,X1))
title(['Rigidezza dinamica per unit� di superficie vs picco della forza'])
xlabel(['Force [N]'])
ylabel(['Apparent Stiffness S [N/m^3]'])
figure
plot (X, mdl.WeightedResiduals.Studentized,'*')
%%
prominanza=1;
distanza=10;

[pks,locs,wi,pi] = findpeaks(temp_zoom, f_zoom,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
    
figure, hold on
plot(f_zoom, temp_zoom)
plot(locs, pks,'*')

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
