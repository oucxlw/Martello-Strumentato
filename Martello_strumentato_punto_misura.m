% Martello strumentato
% Versione per elzborare una singola registrazione
%
set (0,'DefaultFigureWindowStyle','docked')
clc
close all

%%martellatore
clear variables

%%
conf = []; %#ok<NASGU>

%In configuration mettiamo degli indici che ci dicono il campione
%utilizzato, la piastra di carico, la superficie d'appoggio l'adesivo e
%la punta del martello

campione={'massarosa_b'};
piastra={'pesante1'};
appoggio={'nessuno'};
adesivo={'gesso'};
punta={'gomma'};
martellatore=2;
accelerometro=0;
conf = table(campione, piastra, appoggio, adesivo, punta, martellatore, accelerometro) %#ok<NOPTS>

%%
set (0,'DefaultFigureWindowStyle','docked')
clear variables

clc
close all
%clear variables

dati = load ('dati1');
dati(2) = load ('dati2');
dati(3) = load ('dati3');
% dati(4) = load ('dati1a');
% dati(5) = load ('dati2a');
% dati(6) = load ('dati3a');
% dati(7) = load ('dati1b');
% dati(8) = load ('dati2b');
% dati(9) = load ('dati3b');
% dati(10) = load ('dati1c');
% dati(11) = load ('dati2c');
% dati(12) = load ('dati3c');

N = length (dati);
conf.conf = dati(1).conf %#ok<NOPTS>
for i= 2:N
    conf(i).conf = dati(i).conf;
end


fs=dati(1).fs;

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf(1).conf,piastre);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
misure = cell(length(dati),1);
sng = struct('x',misure,'y',[]);
N = length (dati);
for i = 1:N
    sng(i).x = dati(i).data(fs:end-fs,1);
    sng(i).y = dati(i).data(fs:end-fs,conf(i).conf.accelerometro+2);
end

%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri di controllo
%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri fisici
[g,div_F,div_A,fs] = parametri_fisici();
[cal_f, cal_a] = parametri_calibrazione();
% Parametri di ricerca
fine = cell(1,N);
for i = 1:N%length(sng)
    [soglia,delay,inizio,fine{i}] = parametri_ricerca_picchi(fs,sng(i).x);
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

misura = cell2mat(['Campione ',conf(1).conf.campione,', martellatore ',...
    num2str(conf(i).conf.martellatore),', punta ',conf(i).conf.punta,', piastra ',conf(i).conf.piastra,...
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
    sng(i).x = cal_f * div_F * sng(i).x;
    sng(i).y = cal_a * g * div_A * sng(i).y;
end

% Plot dei segnali
for i=1:N
    L = length(sng(i).x);
    dt=1/fs; time=(0:dt:L/fs-dt);
    figure (i)
    subplot(2,1,1), hold on, plot (time, sng(i).x)
    subplot(2,1,2), hold on, plot (time, sng(i).y),
    hold off
end

%<<<<<<<<<<<<<<<<<<<<
% Ricerca dei PICCHI
%<<<<<<<<<<<<<<<<<<<<
prominanza = 10;%25;
distanza = 0.7;
larghezza = 15;
fpass = 35 %#ok<NOPTS>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Sostituzione del segnale con il segnale filtrato
% Attivare solo se � presente troppo rumore

% for i=1:N
%     sng(i).x = x_Hpass{i};
% end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%<<<<<<<<<<<<<<<
% Ricerca Colpi
%>>>>>>>>>>>>>>>

picchi_t = cell(1,N);
picchi1 = cell(1,N);
p = cell(1,N);
w = cell(1,N);
n_picchi = cell(1,N);

filtro = 0 %#ok<NOPTS>
if filtro == 0
    for i=1:N
        [picchi1{i}, n_picchi{i}] = trovacolpi(sng(i).x, soglia, delay, inizio, fine{i});
        n_picchi{i} %#ok<NOPTS>
    end
    
else
    x_Hpass = cell(1,N);
    for i=1:N
        x_Hpass{i} = highpass(sng(i).x,fpass,fs);
    end
    for i=1:N
        
        % [picchi1{i}, n_picchi{i}] = trovacolpi(x, soglia, delay, inizio, fine);
        % n_picchi{i}
        
        [pks1,picchi_t{i},w{i},p{i}] = findpeaks(x_Hpass{i}(1:end),fs,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
        
        picchi_t{i}=picchi_t{i}(w{i}<larghezza);
        
        pks1=pks1(w{i}<larghezza);
        picchi1{i} = picchi_t{i}*fs; % posizioni dei picchi = tempo * fs
        n_picchi{i} = length(picchi1{i}) %#ok<NOPTS>
        
        figure (i)
        subplot(2,1,1), hold on, plot(picchi_t{i},pks1,'*')
        findpeaks(x_Hpass{i},fs,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
    end
end

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
picchi_sel1 %#ok<NOPTS>

[r,~] = size(sng(1).F(:,1));

%<<<<<<<<<<<<<<
% Finestratura
%<<<<<<<<<<<<<<
% Generazione finestra Accelerazione
M_win = round(L_pre/2); % M_win � l'ordine della finestratura, la lunghezza del fronte di salita + discesa
L_plateau = 0.12;%0.0625; % L_plateau indica la lunghezza del plateau rispetto alla massima possibile
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
dt = 1/fs; time1 = 1000*(0:dt:r/fs-dt);
% time1 = (0:dt:r/fs-dt)

% Generazione finestra Forza
L_win_F = L_pre+round(2.25*fs/1000)+M_win/2;
win_F = [ones(L_win_F-(M_win/2),1); curve(end/2:end-1); zeros(Lsample-L_win_F,1)];
win_F(1,1) = 0;

% % Plot delle finestre
% figure (101),hold on
%
% subplot(2,2,2)
% yyaxis right
% plot(time1, win_A )%* 20*log10(max(max(sng(1).A))));
%
% subplot(2,2,1);
% yyaxis right
% plot(time1, win_F )%* max(max(sng(1).F)));
% hold off

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
% figure (15), hold on, plot(10*log10(flip(A_reverse)))

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
win_1 = ones(size(win_F)); % finestra unitaria per non far calcolare la normalizzazione a periodogram

% Periodogram non gestisce pi� di 40 colonne contemporaneamente quindi
% eseguo le operazioni ciclicamente
for i=1:N
    
    [PSD(i).F, ~] = periodogram(sng(i).F_filt, win_1, r, fs); %PSD Forza [N^2]
    [pxy, f] = cpsd(sng(i).F_filt, sng(i).A_filt, win_1, [], r, fs);
    PSD(i).A = abs(pxy).^2;
    %[PSD(i).A, f] = periodogram(sng(i).A_filt, win_1, r, fs); %PSD Accelerazione [g^2]
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
picchi_sel2 = picchi_sel1 - scarti %#ok<NOPTS>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi statistica pre filtraggio
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% bin=[];
% Y=[];
% E=[];
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
picchi_sel2 = picchi_sel2 - scarti %#ok<NOPTS>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo delle grandezze fisiche delle singole misure
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

for i = 1:N
    % psd
    pxx = cpsd (sng(i).F_filt, sng(i).F_filt, win_1, [], r, fs);
    pxx_av = mean (pxx,2);
    PSD(i).Fav = pxx_av;
    %PSD(i).Fav = mean(PSD(i).F, 2);
    
    [pxy, ~] = cpsd(sng(i).F_filt, sng(i).A_filt, win_1, [], r, fs);
    pxy_av = mean(pxy,2);
    PSD(i).Aav = pxy_av.^2;
    %PSD(i).Aav = mean(PSD(i).A, 2);
    
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


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo delle grandezze fisiche mediate sulle N misure
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Calcolo tramite periodogram
pnt_misura.PSD = struct('F',[],'A',[]);
temp_F =[];
temp_A =[];
for i = 1:N
    temp_F = [temp_F, PSD(i).Fav];
    temp_A = [temp_A, PSD(i).Aav];
end

pnt_misura.PSD.F = mean(temp_F,2);
pnt_misura.PSD.A = mean(temp_A,2);

pnt_misura.PSD.V = pnt_misura.PSD.A./(1i*2*pi*f).^2; % velocità
pnt_misura.PSD.D = pnt_misura.PSD.V./(1i*2*pi*f).^2; % displacement

pnt_misura.PSD.K = pnt_misura.PSD.F./pnt_misura.PSD.D; % Dynamic stiffness FFT
pnt_misura.PSD.MI = pnt_misura.PSD.F./pnt_misura.PSD.V; % Mechanical Impedence PSD

% Calcolo tramite fft
pnt_misura.FFT = struct('F',[],'A',[]);

pnt_misura.FFT.F = mean(vertcat([FFT(1:N).F]),2);
pnt_misura.FFT.A = mean(vertcat([FFT(1:N).A]),2);

pnt_misura.FFT.V = pnt_misura.FFT.A./(1i*2*pi*f_fft).^2; % velocità
pnt_misura.FFT.D = pnt_misura.FFT.V./(1i*2*pi*f_fft).^2; % displacement

pnt_misura.FFT.K = pnt_misura.FFT.F./pnt_misura.FFT.D; % Dynamic stiffness FFT
pnt_misura.FFT.MI = pnt_misura.FFT.F./pnt_misura.FFT.V; % Mechanical Impedence PSD


% figure (111), hold on
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


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi delle singole misure
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

m = zeros(N,1);
h = zeros(N,1);
s = zeros(N,1);

for i=1:N
    m(i) = piastre.massa(conf(i).conf.piastra); %massa della piastra in uso
    h(i) = campioni.h(conf(i).conf.campione);
    s(i) = pi*(campioni.d(conf(i).conf.campione)/2)^2;
end

% [R,C]=size(PSD_F);
F_max = [];
fmax_misura = []; % la collezione della frequenza massima della bandwidth della forza a -10dB

% Struttura con i risultati dei singoli colpi
result.indici_colpo.max_A = struct('fd', misure, 'K0', misure,'S_star', [], 'E', []);
result.indici_colpo.min_K = struct('fd', misure, 'K0', misure,'S_star', [], 'E', []);
result.indici_colpo.UNI_10570_1997 = struct('fd', misure, 'nu_d', [], 'K0', [],'S_star', [], 'E', []);
result.indici_colpo.max_V = struct('MI', misure);

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

% Ciclo sulle misure
for mis = 1:N
    
    % Calcolo dimensioni del vettore della Misura
    [RR,CC]=size(PSD(mis).F);
    
    if CC >= 1 % Controllo sulla misura non nulla
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Media della Forza massima impressa nei colpi della misura
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        PSD(mis).devst_F = std (PSD(mis).F, 0, 2);
        PSD(mis).devst_A = std (PSD(mis).A, 0, 2);
        %         PSD(mis).devst_F = std ((PSD(mis).F).^(1/2), 0, 2);
        %         PSD(mis).devst_A = std ((PSD(mis).A).^(1/2), 0, 2);
        
        dF = PSD(mis).devst_F;
        dA = PSD(mis).devst_A;
        b  = abs((1i*2*pi*f).^4);
        FF = abs(PSD(mis).Fav);
        AA = abs(PSD(mis).Aav);
        C  = abs(sng(mis).Cxy(1:1+end/2)).^2;
        KK = PSD(mis).Kav; %#ok<NASGU>
        
        dK = ( ((b./AA).*dF).^2 + ((-b.*FF./AA.^2).*dA).^2 + C.*(b./AA).*(-b.*FF./AA.^2).*dF.*dA  ).^(1/2);
        PSD(mis).devst_K = dK;
        
        
        %PSD(mis).devst_F = std (abs(PSD(mis).F), 0, 2);
        %PSD(mis).devst_A = std (abs(PSD(mis).A), 0, 2);
        
        
        %         dF = PSD(mis).devst_F;
        %         dA = PSD(mis).devst_A;
        %         b  = abs((1i*2*pi*f).^4);
        %         FF = abs(PSD(mis).Fav);
        %         AA = abs(PSD(mis).Aav);
        %         C  = abs(sng(mis).Cxy(1:1+end/2));
        %         KK = PSD(mis).Kav;
        %
        %         dK = ( ((b./AA).*dF).^2 + ((-b.*FF./AA.^2).*dA).^2 + C.*(b./AA).*(-b.*FF./AA.^2).*dF.*dA  ).^(1/2);
        %         PSD(mis).devst_K = dK;
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima della forza media
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 10*log10(PSD(mis).Fav);
        fmax = find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax = f(fmax(1)+f0);
        PSD(mis).fmax = fmax;
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot di segnali e spettri (PSD)
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (101), grid on,
        sgtitle(misura)
        for i = 1:c
            subplot(2,2,1), hold on % Forza nel tempo
            plot (time1, sng(mis).F(:,i),'color',string(colore(mis,:))),
            
            subplot(2,2,2), hold on % Logaritmo dell'accelerazione nel tempo
            % plot (time1, 20*log10(abs(sng(mis).A(:, i))),'color',string(colore(mis,:)))
            plot (time1, sng(mis).A(:, i),'color',string(colore(mis,:)))
            
            subplot(2,2,3), hold on % PSD della Forza
            semilogx (f, 10*log10(abs(PSD(mis).F(:, i))),'--','color',string(colore(mis,:)))
            
            subplot(2,2,4), hold on % PSD dell'accelerazione
            semilogx (f, 10*log10(abs(PSD(mis).A(:, i))),'--','color',string(colore(mis,:)))
        end
        
        % Plot spettri medi
        subplot(2,2,3), hold on % Plot Forza media misura
        plot (f, 10*log10(abs(PSD(mis).Fav)), 'color', string(colore(mis,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        subplot(2,2,4), hold on % Plot Accelerazione media misura
        plot (f, 10*log10(abs(PSD(mis).Aav)), 'color', string(colore(mis,:)), 'LineWidth', 3),
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        % plot frequenza massima sulla PSD della forza
        % subplot (2,2,3), hold on
        % xl=xline(fmax,'.',['F max: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot grafico riassuntivo con dstiff medi di tutti le misure
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (107),hold on,
        
        subplot(4,1,1), hold on % plot coerenza della misura
        plot (f, C,'color',string(colore(mis,:)), 'LineWidth', 2,...
            'DisplayName',['Misura N�',num2str(mis)]),
        
        subplot(4,1,[2,3]), hold on % Plot rigidezza dinamica della misura
        plot (f, 10*log10(PSD(mis).Kav),'color',string(colore(mis,:)),...
            'LineWidth', 2, 'DisplayName',['Misura N�',num2str(mis)]),
        
        % Plot deviazione standard sulla misura
        %semilogx (f, 10 * log10( PSD(mis).Kav ) + 10/log(10) * dK./PSD(mis).Kav, '-','color',string(colore(mis,:)), 'LineWidth', 1),
        %semilogx (f, 10 * log10( PSD(mis).Kav ) - 10/log(10) * dK./PSD(mis).Kav, '-','color',string(colore(mis,:)), 'LineWidth', 1),
        
        plot (f, 10*log10(PSD(mis).Kav)+10*log10(1+dK./PSD(mis).Kav),...
            '-.','color',string(colore(mis,:)),'LineWidth',1,'HandleVisibility','off'),
        plot (f, 10*log10(PSD(mis).Kav)-10*log10(1+dK./PSD(mis).Kav),...
            '-.','color',string(colore(mis,:)),'LineWidth',1,'HandleVisibility','off'),
        
        % semilogx (f, 10 * log10( PSD(mis).Kav + dK ), '-.','color',string(colore(mis,:)), 'LineWidth', 1),
        % semilogx (f, 10 * real( log10( PSD(mis).Kav - dK ) ), '-.','color',string(colore(mis,:)), 'LineWidth', 1),
        
        % Plot limite in frequenza
        xline(PSD(mis).fmax, '.', ['Limite in frequenza: ', ...
            num2str(round(PSD(mis).fmax)),' Hz'],'color',string(colore(mis,:)),'HandleVisibility','off');
        %Plot fase
        subplot(4,1,4), hold on % Plot fase
        semilogx (f_fft, 180/pi * angle(FFT(mis).Kav), 'color',...
            string(colore(mis,:)),'LineWidth', 2, 'DisplayName',['Misura N�',num2str(mis)])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot Accelerazione e forza
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        figure (110), hold on
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
        
        lim_sup =  find(f > 500); %min(vertcat([PSD.fmax]))*1);%0.35);
        lim_inf =  find(f < 10);
        
        % secondo minimo di K
        min_k = min(PSD(mis).Kav(lim_inf(end):lim_sup(1)));
        fd = f( lim_inf(end)+find(PSD(mis).Kav(lim_inf(end):lim_sup(1)) == min_k) );
        result.indici.min_K(mis).fd = fd;
        result.indici.min_K(mis).K0 = (2*pi*fd)^2 * m(mis);
        result.indici.min_K(mis).S_star = (2*pi*fd)^2 * m(mis)/s(mis);
        result.indici.min_K(mis).E  = (2*pi*fd)^2 * m(mis) * h(mis)/s(mis);
        
        % secondo massimo di A
        temp = PSD(mis).Aav;%./PSD(mis).Fav;
        max_a = max(abs(temp(lim_inf(end):lim_sup(1))));
        fd = f(lim_inf(end)+find (abs (temp(lim_inf(end):lim_sup(1))) == max_a));
        result.indici.max_A(mis).fd = fd;
        result.indici.max_A(mis).K0 = (2*pi*fd)^2 * m(mis);
        result.indici.max_A(mis).S_star = (2*pi*fd)^2 * m(mis)/s(mis);
        result.indici.max_A(mis).E  = (2*pi*fd)^2 * m(mis) * h(mis)/s(mis);
    end
end
clear x_Hpass;

%<<<<<<<<<<<<<<<<<<<<<<<<<
%Salvataggio degli output
%<<<<<<<<<<<<<<<<<<<<<<<<<
frequenze.f = f;
frequenze.f_fft = f_fft;

% Figura Segnali e Spettri

% Plot delle finestre
figure (101),hold on

subplot(2,2,2)
yyaxis right
plot(time1, win_A,'LineWidth',2)%* 20*log10(max(max(sng(1).A))));
ylabel('Finestra')
ylim([-0.1, 1.1]), yticks (0:0.1:1)

subplot(2,2,1);
yyaxis right
plot(time1, win_F,'LineWidth',2)%* max(max(sng(1).F)));
ylabel('Finestra')
hold off
ylim([-0.1, 1.1]), yticks (0:0.1:1)

% Settaggio parametri grafico
figure (101), hold on
subplot(2,2,1), hold on
yyaxis left
xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Forza')
grid on, %xlim([0 10])
subplot(2,2,2), hold on,
yyaxis left
xlabel('Time [ms]'), ylabel('Amplitude 20 log10 ([m/s^2])'), title('Accelerazione'),
grid on, %xlim([0 10])
subplot(2,2,3), hold on
xlabel('log(Frequency) [Hz]'), ylabel('10 log |PSD| (dB ref 1 N/Hz)'), title('PSD Forza')
grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
subplot(2,2,4), hold on,
xlabel('log(Frequency) [Hz]'), ylabel('10 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Accelerazione'),
grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
% Salvataggio
figure (101),saveas (gcf, 'Segnali e spettri.fig')

% Trazformata calcolata in un intervallo di frequenze ridotto
f1=100 %#ok<NOPTS>
f2=1000 %#ok<NOPTS>
% f_zoom = f1:0.5:f2;
%
% for i=1:N
%     %[PSD_Ftemp, ~, pxxc] = periodogram( sng(i).F_filt, win_1, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Forza [N^2]
%     %[PSD_Atemp, f_2, pyyc] = periodogram(sng(i).A, win_1, f_zoom, fs, 'ConfidenceLevel',0.60); %PSD Accelerazione [g^2]
%
%     [pxx, ~] = cpsd(sng(1).F_filt, sng(i).F_filt, win_1, [], f_zoom, fs);
%     [pxy, ~] = cpsd(sng(1).F_filt, sng(i).A_filt, win_1, [], f_zoom, fs);
%
%     PSD(i).F_zoom = mean (pxx, 2);
%     PSD(i).A_zoom = mean (pxy, 2).^2;
%
%     % PSD(i).F_zoom = mean(PSD_Ftemp, 2);
%     % PSD(i).A_zoom = mean(PSD_Atemp, 2);
% end

% Stampa a schermo di uno zoom della trasformata ma sembra inutile

% figure (107),
% subplot(4,1,[2,3]), hold on
% for i=1:N
%     plot (f_zoom,10*log10(abs(PSD(i).F_zoom./PSD(i).A_zoom.*(1i*2*pi*f_zoom').^4)),'LineWidth', 1 ,'color',string(colore(i,:)) )
% end
% % sgtitle(['PSD della rigidezza dinamica [N/mHz]. Campione:',cell2mat([conf(i).conf.campione,...
% %     ' + ',conf(i).conf.piastra,'. Adesivo: ',conf(i).conf.adesivo,'.'])])
% grid on, hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot accelerazione normalizzata e Zoom
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%Accelerazione normalizzata
% figure (400), hold on
% for i=1:N
%     plot(f,10*log10(PSD(i).Aav./PSD(i).Fav), 'color',string(colore(i,:)))
%
%     % psd zoom
%     plot(f_zoom, 10*log10(PSD(i).A_zoom./PSD(i).F_zoom), '+', 'color',string(colore(i,:)))
% end
%
% %plot(f_zoom, 10*log10(pyyc./PSD(i).F_zoom), '-*')
% grid on
% set(gca, 'XScale', 'log'),
% xlim([10 2000])
% title('Accelerazione normalizzata','FontSize',18)
% xlabel('Frequency [Hz]','FontSize',12),
% ylabel('Acceleration/Force [Kg^{-1}]','FontSize',12),
% saveas(gcf,'Accelerazione.fig');

% accelerazione/forza
% figure (401), hold on
% for i=1:3
%     plot(f,PSD(i).Aav ./ PSD(i).Fav)
% end
% grid on
% set(gca, 'XScale', 'log')
% xlim([10 1000])


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Figura confronto Accelerazione-Coerenza secondo Vasquez
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


figure (405), hold on
title ('Frequenza di risonanza (V.F. Vazquez, S.E. Paje, 2012) ','FontSize',16)
f_analisi=[20 1000];
yyaxis left
hold on
massimi = cell(1,N);
%fr_misura = zeros(1,N);
in = find(f>f_analisi(1),1);
out = find(f<f_analisi(2));
out = out(end);

for i = 1:N
    temp = abs(PSD(i).A);
    temp_massimi=zeros(1,size(temp,2));
    for j=1:size(temp,2)
        temp_massimi(j) = f(temp(:,j)==max(temp(in:out,j)));%cerco la frequanza di risonanza
    end
    %massimi(i) = f(temp==max(temp(in:out)));%cerco la frequanza di risonanza
    result.indici_colpo.max_V(i).fd = temp_massimi;
    %pnt_misura.indici.max_A(i).fd = massimi(i); %scrivo il risultato
    plot(f,100*(temp/max(temp(in:out))),'DisplayName',['Sequenza n�',...
        num2str(i)],'LineWidth',2)
end
ylim ([0 150])
ylabel('Accelerazione [%]','FontSize',18),

yyaxis  right
hold on
for i = 1:N
    plot(f_fft,abs(sng(i).Cxy),'DisplayName',['Sequenza n�',num2str(i)],...
        'LineWidth',2)
end
ylabel('Coerenza','FontSize',18),
xlim (f_analisi),ylim([0 1]);
S_primo_0 = zeros(1,N);
figure (406) %S' in funzone della forza
for i=1:N
    result.indici_colpo.max_V(i).K0 = (2*pi* result.indici_colpo.max_V(i).fd ).^2 * m(i);
    result.indici_colpo.max_V(i).S_star = (2*pi* result.indici_colpo.max_V(i).fd ).^2 * m(i)/s(i);
    result.indici_colpo.max_V(i).E  = (2*pi* result.indici_colpo.max_V(i).fd ).^2 * m(i) * h(i)/s(i);
    plot(max(sng(i).F),  result.indici_colpo.max_V(i).S_star/10^6,'+','LineWidth',2,...
        'color', string(colore(i,:)),'DisplayName', ['Sequenza n�', num2str(i)])
    
    X = max(sng(i).F);
    Y = result.indici_colpo.max_V(i).S_star;
    P = polyfit(X,Y,1);
    
    X_exp = 1:max(X*1.5);
    Y_fit = P(1)*X_exp+P(2);
    S_primo_0(i) = Y_fit(1);
    hold on;
    plot(X_exp,Y_fit/10^6,'r-.', 'color', string(colore(i,:)),...
        'DisplayName',['Fit lineare n�', num2str(i)]);
end
grid on;
ylabel('Digidezza Dinamica apparente x10^6 [N/m]','FontSize',18),
xlabel ('Forza [N]', 'FontSize', 18)
title ('Rigidezza dinamica in funzione della forza impressa','FontSize',20)
save('S_primo_0','S_primo_0');
exportgraphics(gcf,'Dstiff_vs_F-Vasquez.pdf', 'BackgroundColor', 'none')
saveas(gcf, cell2mat(['Dstiff_vs_F-Vasquez_',conf(i).conf.campione,'_',conf(i).conf.piastra,'.fig']))


figure (402), hold on
title ('Frequenza di risonanza (V.F. Vazquez, S.E. Paje, 2012) ','FontSize',16)

yyaxis left
hold on
fr_misura = zeros(1,N);
for i = 1:N
    temp = abs(PSD(i).Aav);
    in = find(f>f_analisi(1),1);
    out = find(f<f_analisi(2));
    out = out(end);
    fr_misura(i) = f(temp==max(temp(in:out)));%cerco la frequanza di risonanza
    pnt_misura.indici.max_A(i).fd=fr_misura(i); %scrivo il risultato
    plot(f,100*(temp/max(temp(in:out))),'DisplayName',['Sequenza n�',...
        num2str(i),' f_r=',num2str(fr_misura(i)),' Hz'],'LineWidth',2)
end
ylim ([0 150])
ylabel('Accelerazione [%]','FontSize',18),

yyaxis  right
hold on
for i = 1:N
    plot(f_fft,abs(sng(i).Cxy),'DisplayName',['Sequenza n�',num2str(i)],...
        'LineWidth',2)
end
ylabel('Coerenza','FontSize',18),

xlim (f_analisi),ylim([0 1]);
xlabel('Frequenza [Hz]','FontSize',18),
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
legend('FontSize',12)
% Requires R2020a or later
exportgraphics(gcf,'Acc-Coher-Vasquez.pdf', 'BackgroundColor', 'none')
saveas(gcf, cell2mat(['Acc-Coher-Vasquez_',conf(i).conf.campione,'_',conf(i).conf.piastra,'.fig']))

for i=1:N
    %pnt_misura.indici.max_A(i).K0 = fr_misura(i);
    pnt_misura.indici.max_A(i).K0 = (2*pi*fr_misura(i))^2 * m(i);
    pnt_misura.indici.max_A(i).S_star = (2*pi*fr_misura(i))^2 * m(i)/s(i);
    pnt_misura.indici.max_A(i).E  = (2*pi*fr_misura(i))^2 * m(i) * h(i)/s(i);
end
s_pnt_misura = vertcat(pnt_misura.indici.max_A.S_star);
save('S_primo','s_pnt_misura');



%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dstiff totale e settaggio parametri grafico
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (107)
sgtitle(['PSD della rigidezza dinamica [N/mHz]. Campione:',cell2mat([conf(i).conf.campione,...
    ' + ',conf(i).conf.piastra,'. Adesivo: ',conf(i).conf.adesivo,'.'])])
% sgtitle('PSD della rigidezza dinamica [N/mHz] media su tutte le misure.')

F_max_av = zeros(N,1);
for i=1:N
    F_max_av(i) = max(max(sng(i).F));
end

subplot(4,1,1), hold on
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),ylim([0 1])
xlabel('Frequenza [Hz]'), ylabel('Coerenza'),
% legend(['misura 1 (',num2str(round(F_max_av (1))),' N)'],...
%     ['misura 2 (',num2str(round(F_max_av (2))),' N)'],...
%     ['misura 3 (',num2str(round(F_max_av (3))),' N)'])
legend ('FontSize', 12)

% Plot della massa della piastra
subplot(4,1,[2,3]), hold on
% Plotto anche le masse delle altre misure se hanno valori diversi da m(1)
for mis=1:N
    if mis == 1
        plot(f,10*log10(m(mis)*(2*pi*f).^4),'--k','DisplayName','Massa del campione');
    else
        if m(mis)~=m(1)
            plot(f,10*log10(m(mis)*(2*pi*f).^4),'--k','DisplayName',...
                'Massa del N�',num2str(i),' campione');
        end
    end
end
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),ylim([80 200])
xlabel('Frequenza [Hz]'), ylabel('PSD Rigidezza Dinamica [dB @ 1 N/m]'),
legend ('FontSize', 12)

subplot(4,1,4), hold on
set(gca, 'XScale', 'log'), yticks([-180 -90 0 90 180])
grid on, xlim([ascissamin ascissamax]),ylim([-180 180])
xlabel('Frequenza [Hz]'), ylabel('Fase [ang]'),
legend ('FontSize', 12)

% legend(['misura 1 (',num2str(round(F_max_av (1))),' N)'],...
%     ['misura 2 (',num2str(round(F_max_av (2))),' N)'],...
%     ['misura 3 (',num2str(round(F_max_av (3))),' N)'])

% Salvataggio
saveas(gcf, cell2mat(['Collezione Dstiff_',conf(i).conf.campione,'_',conf(i).conf.piastra,'bw',num2str(bandwidth),'.fig']))
exportgraphics(gcf, cell2mat(['Collezione Dstiff_',conf(i).conf.campione,'_',...
    conf(i).conf.piastra,'bw',num2str(bandwidth),'.pdf']), 'BackgroundColor', 'none')

%saveas(gcf, cell2mat(['Collezione Dstiff_',conf.campione,'_',conf.piastra,'.png']))

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Figura Accelerazione e Forza
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (110), hold on

subplot(2,1,1),hold on
set(gca, 'XScale', 'log')%, set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Forza [dB @ 1 N]'),

subplot(2,1,2),hold on,
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
grid on, xlim([ascissamin ascissamax]),%ylim([120 220])
xlabel('Frequenza [Hz]'), ylabel('PSD Accelerazione [dB @ 1 m/s^2]'),
saveas(gcf, cell2mat(['Forza_vs_Accelerazione_',conf(i).conf.campione,'_',conf(i).conf.piastra,'.fig']))

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo delle grandezze sintetiche della misura
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%F_filtall = reshape(F_filt, [],1); % Qui non si capisce bene che senso abbia
%A_filtall = reshape(A_filt, [],1);
%[r,c]=size(F_filt);

% Calcolo coerenza
%[Cxy_sintetico,f_coer] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],f_fft,fs); % questa coerenza non sembra servire a qualcosa
% forse vale la pena di usare la coerenza media

Cxy_sintetico = mean (vertcat([sng(1:N).Cxy]), 2);



% deviazione standard
pnt_misura.PSD.devst_K = mean(vertcat([PSD(1:N).devst_K]).^2, 2);
pnt_misura.PSD.devst_F = std (vertcat([PSD(1:N).F]),0,2);
pnt_misura.PSD.devst_A = std (vertcat([PSD(1:N).A]),0,2);

dF = pnt_misura.PSD.devst_F;
dA = pnt_misura.PSD.devst_A;
b  = abs((1i*2*pi*f).^4);
FF = pnt_misura.PSD.F;
AA = pnt_misura.PSD.A;
C  = abs(Cxy_sintetico(1:1+end/2));
KK = pnt_misura.PSD.K;
% K = F*(1i*2*pi).^4 /A

dK = ( ((b./AA).*dF).^2 + ((-b.*FF./AA.^2).*dA).^2 + C.*(b./AA).*(-b.*FF./AA.^2).*dF.*dA  ).^(1/2);

pnt_misura.PSD.devst_K_sintetico = dK;
% Salvataggio deviazione standard sintetica rigidezza dinamica PSD
save ('DevSt_sintetica', 'dK');
save ('DevSt_sintetica', 'dK');


% % Salvataggio delle k0 misurate bin per bin
% save (cell2mat(['Dstiffness_bin_',conf.campione,'_',conf.piastra,'_Fft.mat']),'PSD_K_bin');


figure (106)
grid on %xlim([0 10])
set(gca, 'XScale', 'log')
xlim([10 20000])
hold on

for i=1:3
    subplot(2,1,1),hold on
    title('Impedenza meccanica')
    plot(f,10*log10(abs(PSD(i).Fav ./PSD(i).Vav)));
    subplot(2,1,2),hold on
    plot(f_fft,180/pi*unwrap(angle(FFT(i).Fav./FFT(i).Vav)))
    
end
subplot(2,1,1),hold on, set(gca, 'XScale', 'log'), grid on

subplot(2,1,2),hold on, set(gca, 'XScale', 'log'), grid on
ylim([-200 200])

% % velocit�/forza
% figure , hold on
% plot(f,PSD_Vav_misura(:,1)./PSD_Fav_misura(:,1))
% plot(f,PSD_Vav_misura(:,2)./PSD_Fav_misura(:,2))
% plot(f,PSD_Vav_misura(:,3)./PSD_Fav_misura(:,3))
% grid on
% set(gca, 'XScale', 'log')
% xlim([10 2000])
%
% % spostamento/forza
% figure , hold on
% plot(f,PSD_Dav_misura(:,1)./PSD_Fav_misura(:,1))
% plot(f,PSD_Dav_misura(:,2)./PSD_Fav_misura(:,2))
% plot(f,PSD_Dav_misura(:,3)./PSD_Fav_misura(:,3))
% grid on
% set(gca, 'XScale', 'log')
% xlim([10 2000])

%
% %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% % Calcolo frequenza di risonanza e K
% %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% m=piastre.massa(conf(i).conf.piastra); %massa della piastra in uso
% h=campioni.h(conf(i).conf.campione);
% s=pi*(campioni.d(conf(i).conf.campione)/2)^2;
% lim_sup =  find(f > min(vertcat([PSD.fmax]))*0.7);
% lim_inf =  find(f < 100);
% min_k = min(PSD_Kav_pgram(lim_inf(end):lim_sup(1)));
%
% % basato su minimo Dynamic Stiffness
% fd = f( lim_inf(end)+find(PSD_Kav_pgram(lim_inf(end):lim_sup(1)) == min_k)); % calcolata con la fft meno
% K0_av = (2*pi*fd)^2*m %#ok<NOPTS>
% S_star_av = K0_av/s %#ok<NOPTS>
% E_av = K0_av*h/s %#ok<NOPTS>
%
% %basato sul massimo dell'accelerazione/forza
% fr_A = zeros(N,1);
% for i=1:N
%     temp = PSD(i).Aav./PSD(i).Fav;
%     max_A = max(temp(lim_inf(end):lim_sup(1)));
%     fr_A(i) = f(lim_inf(end)+find(temp(lim_inf(end):lim_sup(1)) == max_A)); % calcolata con la fft meno
% end
% K0_av_A = (2*pi*mean(fr_A))^2*m
% S_star_av_A = K0_av_A/s
% E_av_A = K0_av_A*h/s
%
% figure (5), hold on
% xl=xline (f(lim_sup(1)),'.',('limite superiore')); xl.LabelVerticalAlignment = 'bottom';
% xl=xline (f(lim_inf(end)),'.',('limite superiore')); xl.LabelVerticalAlignment = 'bottom';
%
% %figure (107),hold on, plot(f,20*log10(abs(K(1:length(f)))),'r.','LineWidth', 2)
%
% %<<<<<<<<<<<<<<<<<<<<<<<<<
% % Salvataggio dati misura
% %<<<<<<<<<<<<<<<<<<<<<<<<<
%
% save (['misura'],'Cxy','PSD_Kav_misura','FFT_K_misura','devst_K_misura','devst_K_sintetico',...
%     'S_av_bin','conf','piastre','campioni','S_star_av_A')
%
% save('frequenze.mat','f','f_fft')
%
% save (cell2mat(['Dstiffness_bin_',conf(i).conf.campione,'_',conf(i).conf.piastra,'_Fft.mat']),'PSD_K_bin');

save('Outputs.mat', 'sng', 'PSD', 'FFT', 'result','frequenze', 'win_1', 'win_A', 'conf');
%%
%<<<<<<<<<<<<<<<<<<<<<<
% Mechanical impedence
%<<<<<<<<<<<<<<<<<<<<<<

% Calcolo dell'impedenza meccanica come suggerito da
% Pavement stiffness measurements in relation to mechanical impedance
% Mingliang Li, Wim van Keulen, Halil Ceylan, Dongwei Cao,
% mechanical_impedence plotta anche grafici di velocit� vs accelerazione
% e velocit� vs forza

[tempstruct,mdls,MI_av,MI_av_dev,MI_lin,MI_lin_R,MI_exp,MI_exp_R] = mechanical_impedence(sng, N,fs);
MI_mean = zeros(1,N);
% Applico a result le MI calcolate in mechanical_impedence
for i = 1:N
    result.indici_colpo.max_V(i).MI = tempstruct(i).MI;
    MI_mean(i) = mean(result.indici_colpo.max_V(i).MI);
end
clear tempstruct

for i=1:N
    result.indici_colpo.max_V(i).MI_av = MI_av(i); %mean(result.indici_colpo.max_V(i).MI);
    result.indici_colpo.max_V(i).MI_lin = MI_lin(i);
    result.indici_colpo.max_V(i).MI_exp = MI_exp(i);
    
end
save('MI_modelli_fit.mat','MI_mdls');
save('MI_mean.mat','MI_av');
save('MI_exp.mat','MI_exp');
save('MI_lin.mat','MI_lin');

save('Outputs.mat', 'sng', 'PSD', 'FFT', 'result','frequenze', 'win_1', 'win_A', 'conf');

%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Impedenza meccanica stand alone
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
frequenze = Inputs.frequenze;

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf(1).conf,piastre);

N = length(result.indici_colpo.max_A);
fs = 52100;


for i = 1:N
    m(i) = piastre.massa(conf(i).conf.piastra); %massa della piastra in uso
    h(i) = campioni.h(conf(i).conf.campione);
    s(i) = pi*(campioni.d(conf(i).conf.campione)/2)^2;
end

clear Inputs;

% Plotting
ascissamin = 100;          % Frequenza minima da plottare nei grafici
ascissamax = 2000;       % Frequenza massima da plottare nei grafici

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

%<<<<<<<<<<<<<<<<<<<<<<
% Mechanical impedence
%<<<<<<<<<<<<<<<<<<<<<<

% Calcolo dell'impedenza meccanica come suggerito da
% Pavement stiffness measurements in relation to mechanical impedance
% Mingliang Li, Wim van Keulen, Halil Ceylan, Dongwei Cao,
% mechanical_impedence plotta anche grafici di velocit� vs accelerazione
% e velocit� vs forza

[tempstruct,MI_zero] = mechanical_impedence(sng, N,fs);
MI_mean = zeros(1,N);
% Applico a result le MI calcolate in mechanical_impedence
for i = 1:N
    result.indici_colpo.max_V(i).MI = tempstruct(i).MI;
    MI_mean(i) = mean(result.indici_colpo.max_V(i).MI);
end
clear tempstruct

MI_av = zeros(1,N);
for i=1:N
    result.indici_colpo.max_V(i).MI_av = MI_zero; %mean(result.indici_colpo.max_V(i).MI);
end
save('MI_mean.mat','MI_mean');
save('MI_zero.mat','MI_zero');
save('Outputs.mat', 'sng', 'PSD', 'FFT', 'result','frequenze', 'win_1', 'win_A', 'conf');



%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Analisi risultati singoli colpi
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
[campioni] = tabella_campioni (conf(1).conf,piastre);

N = length(result.indici_colpo.max_A);
fs = 52100;


for i = 1:N
    m(i) = piastre.massa(conf(i).conf.piastra); %massa della piastra in uso
    h(i) = campioni.h(conf(i).conf.campione);
    s(i) = pi*(campioni.d(conf(i).conf.campione)/2)^2;
end

clear Inputs;

% Plotting
ascissamin = 100;          % Frequenza minima da plottare nei grafici
ascissamax = 2000;       % Frequenza massima da plottare nei grafici

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

% accelerazione/forza
lim_sup =  find(f > min(vertcat([PSD.fmax]))*1);
lim_inf =  find(f < 10);
f1 = f(lim_inf(end));

f2 = f(lim_sup(1));
f_zoom = f1:0.2:f2;
prominanza=0;%0.2;
distanza=20;
result.indici_colpo.max_A(1).locs = {1};
result.indici_colpo.max_A(1).A_norm_zoom = [];

for i = 1:N
    [~,cc] = size(PSD(i).F);
    [PSD(i).F_zoom, ~, pxxc] = periodogram(sng(i).F_filt, win_1, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Forza [N^2]
    [PSD(i).A_zoom, ~, pyyc] = periodogram(sng(i).A, win_A, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Accelerazione [g^2]
    
    figure (120), hold on
    for j = 1:cc
        % secondo minimo di K
        PSD(i).K(:,j) = PSD(i).F(:,j) ./ PSD(i).A(:,j) .* (2*pi*f).^4;
        plot (f, 10*log10(PSD(i).K(:,j)), 'color', string(colore(i,:)));
        
        PSD(i).K_zoom(:,j) = PSD(i).F_zoom(:,j) ./ PSD(i).A_zoom(:,j) .* (1i*2*pi*f_zoom').^4;
        plot (f_zoom, 10*log10(PSD(i).K_zoom(:,j)), 'color', string(colore(i+1,:)));
        
        min_k = min(PSD(i).K_zoom(lim_inf(end):lim_sup(1), j));
        fd = f_zoom(lim_inf(end) + find(PSD(i).K_zoom(lim_inf(end):lim_sup(1), j) == min_k) );
        result.indici_colpo.min_K(i).fd(j) = fd;
        result.indici_colpo.min_K(i).K0(j) = (2*pi*fd)^2 * m(i);
        result.indici_colpo.min_K(i).S_star(j) = (2*pi*fd)^2 * m(i)/s(i);
        result.indici_colpo.min_K(i).E(j)  = (2*pi*fd)^2 * m(i) * h(i)/s(i);
        %         plot (f, 10*log10(1.487*1.05*(2*pi*f).^4),'b--', 'LineWidth', 2)
        %         plot (f, 10*log10(1.487*(2*pi*f).^4),'b--', 'LineWidth', 2)
        %         plot (f, 10*log10(1.487*0.95*(2*pi*f).^4),'b--', 'LineWidth', 2)
        plot (f, 10*log10((2*piastre.massa(conf(i).conf.piastra)+0.0573)^2*1.05^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
        plot (f, 10*log10((2*piastre.massa(conf(i).conf.piastra)+0.0573)^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
        plot (f, 10*log10((2*piastre.massa(conf(i).conf.piastra)+0.0573)^2*0.95^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
        plot (f, 10*log10((piastre.massa(conf(i).conf.piastra)+0.0573)^2*1.05^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
        plot (f, 10*log10((piastre.massa(conf(i).conf.piastra)+0.0573)^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
        plot (f, 10*log10((piastre.massa(conf(i).conf.piastra)+0.0573)^2*0.95^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
        
        
        
    end
    xl=xline(f(lim_inf(end)),'.',['Lim inf ',f(lim_inf(end)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    xl=xline(f(lim_sup(1)),'.',['Lim sup ',f(lim_sup(1)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    grid on
    set(gca, 'XScale', 'log')
    xlim([ascissamin ascissamax])
    
    % secondo massimo di A
    figure (130+i), hold on
    for j = 1:cc
        % seconfo massimo di A
        temp = PSD(i).A(:,j) ./ PSD(i).F(:,j);
        plot (f, temp, 'color', string(colore(i,:)));
        
        temp_zoom = PSD(i).A_zoom(:,j) ./ PSD(i).F_zoom(:,j);
        result.indici_colpo.max_A(i).A_norm_zoom(:,j) = temp_zoom;
        
        plot (f_zoom, temp_zoom, 'color', string(colore(i+1,:)));
        
        [pks, locs, wi, pr] = findpeaks(temp_zoom, f_zoom,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents'); %#ok<ASGLU>
        result.indici_colpo.max_A(i).locs{:,j} = locs;
        plot(result.indici_colpo.max_A(i).locs{j}, pks,'*')
        
        max_A = max(temp_zoom);
        if ~isempty(result.indici_colpo.max_A(i).locs{j})
            fd = result.indici_colpo.max_A(i).locs{j}(1);
        end
        if length(result.indici_colpo.max_A(i).locs{j}) > 1
            f2=result.indici_colpo.max_A(i).locs{j}(2);
        end
        
        result.indici_colpo.max_A(i).fd(j) = fd;
        result.indici_colpo.max_A(i).f2(j) = f2;
        result.indici_colpo.max_A(i).K0(j) = (2*pi*fd)^2 * m(i);
        result.indici_colpo.max_A(i).S_star(j) = (2*pi*fd)^2 * m(i)/s(i);
        result.indici_colpo.max_A(i).E(j)  = (2*pi*fd)^2 * m(i) * h(i)/s(i);
    end
    xl=xline(f(lim_inf(end)),'.',['Lim inf ',f(lim_inf(end)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    xl=xline(f(lim_sup(1)),'.',['Lim sup ',f(lim_sup(1)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    grid on
    set(gca, 'XScale', 'log')
    xlim([ascissamin f(lim_sup(1))]), %ylim([0 2*max_A])
    
end

%<<<<<<<<<<<<<<<<<<<<<
% Plot Massa dinamica
%<<<<<<<<<<<<<<<<<<<<<

figure, hold on
for i=1:N
    plot (f ,(PSD(i).F./PSD(i).A).^(1/2),'color', string(colore(i,:)))
end
ylim([0 5])
set(gca, 'XScale', 'log')
xlim([0 2000])
hold on, plot (f ,ones(size(f))*piastre.massa('pesante1') ,'b--', 'LineWidth', 2)
hold on, plot (f ,ones(size(f))*piastre.massa('pesante1')*1.05 ,'b', 'LineWidth', 1)
hold on, plot (f ,ones(size(f))*piastre.massa('pesante1')*0.95 ,'b', 'LineWidth', 1)

hold on, plot (f ,ones(size(f))*piastre.massa('pesante2'),'b--', 'LineWidth', 2)
hold on, plot (f ,ones(size(f))*piastre.massa('pesante2')*1.05 ,'b', 'LineWidth', 1)
hold on, plot (f ,ones(size(f))*piastre.massa('pesante2')*0.95 ,'b', 'LineWidth', 1)

savefig('Calibrazione.fig');

%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Fit bassato sulla ricerca di massimi e minimi
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure (120), hold on,
for i = 1:N
    [~,cc] = size(PSD(i).F);
    for j = 1:cc
        %plot (max(sng(i).F_filt(:,j)), result.indici_colpo.min_K(i).S_star(j)/10^6, '*', 'color', string(colore(i,:)))
        plot (max(sng(i).F_filt(:,j)), result.indici_colpo.max_A(i).S_star(j)/10^6, 'o', 'color', string(colore(i+1,:)), 'lineWidth', 3)
    end
end
grid on, %ylim([140 200])
X = max(vertcat([ sng((1:N)).F_filt]));
Y = vertcat([result.indici_colpo.max_A((1:N)).S_star]) / 10^6;

modelfun = @(b,x)b(1)*exp(-(x-b(4))/b(2))+b(3)
mdl = fitnlm(X, Y, modelfun, [100 50 10 0]); % Modello dei primi picchi
[fit,R] = nlinfit(X, Y, modelfun, [200 150 150 0]);

X_exp = 1:max(X*1.5);
plot (X_exp, feval(mdl,X_exp))
title(('Rigidezza dinamica per unit� di superficie vs picco della forza'))
xlabel(('Force [N]'))
ylabel(('Apparent Stiffness S [N/m^3]'))
figure
plot (X, mdl.WeightedResiduals.Studentized,'*')

k=1;
for i = 1:N
    [~,cc] = size(PSD(i).F);
    for j = 1 : cc
        if length(result.indici_colpo.max_A(i).locs{j}) > 1
            Y2(k) = (2*pi*result.indici_colpo.max_A(i).locs{j}(2))^2*m/s/10^6; %#ok<SAGROW>
            result.indici_colpo.max_A(i).locs{j}(2)
            X2(k) = max(sng(i).F_filt(:,j)); % forza
            k = k+1;
        end
    end
end
figure (120)
hold on
plot (X2, Y2, '*', 'color', string(colore(i,:)))
mdl2 = fitnlm(X2, Y2, modelfun, [100 50 10 0]); % Modello dei secondi picchi
[fit2,R2] = nlinfit(X2, Y2, modelfun, [200 150 150 0]);


X2_exp = 1:max(X2*1.5);
plot (X2_exp, feval(mdl2,X2_exp))
title(('Rigidezza dinamica per unit� di superficie vs picco della forza'))
xlabel(('Force [N]'))
ylabel(('Apparent Stiffness S [N/m^3]'))
figure
plot (X2, mdl2.WeightedResiduals.Studentized,'*')

prominanza=1;
distanza=10;

[pks,locs,wi,p3] = findpeaks(temp_zoom, f_zoom,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');

% figure, hold on
% plot(f_zoom, temp_zoom)
% plot(locs, pks,'*')
%% Plot 3d
%<<<<<<<<<<<<<<
% intestazione
%<<<<<<<<<<<<<<
set (0,'DefaultFigureWindowStyle','docked')
clear variables
close all
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
frequenze = Inputs.frequenze;

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf(1).conf,piastre);

N = length(result.indici_colpo.max_A);
fs = 52100;

for i = 1:N
    m(i) = piastre.massa(conf(i).conf.piastra); %massa della piastra in uso
    h(i) = campioni.h(conf(i).conf.campione);
    s(i) = pi*(campioni.d(conf(i).conf.campione)/2)^2;
end

clear Inputs;

% Plotting
ascissamin = 100;          % Frequenza minima da plottare nei grafici
ascissamax = 2000;       % Frequenza massima da plottare nei grafici

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

%<<<<<<<<<
% Plot 3d
%<<<<<<<<<

%figure (201)

y_max = 500;
y_min = 0;
z_lim = 1;

for i = 1:N
    
    % utilizzo il massimo della forza
    X3 = max(vertcat([sng(i).F_filt]));
    
    % utilizzo dell'energia al posto della forza massima
    %X3 = cumtrapz(    vertcat([sng(i).F_filt.^2]));
    %X3 = X3(end,:);
    
    % utilizzo la Forza RMS
    temp2 = PSD(i).F;
    F_RMS = (sum(temp2,1)/2).^(1/2);
    X3 = F_RMS;
    [X3, Fsort] = sort(X3);
    Y3 = (2*pi*f).^2*m(i)/s(i)/10^6;
    %Y3 = (2*pi*f_zoom).^2*m(i)/s(i);
    %Y3 = f_zoom;
        
    %Z3 = vertcat([ result.indici_colpo.max_A(i).A_norm_zoom]);
    Z3 = vertcat( PSD(i).A ./ max(PSD(i).A((Y3<y_max & Y3>y_min),:)) );
    Z3=abs(Z3);
    Z3 = Z3(:,Fsort);
    % mesh(X3, Y3, Z3)
    % %ylim([0 6*10^8])
    
    figure (i)
    S = mesh(X3,Y3,Z3);
    % Create zlabel
    zlabel('Acc/Force','FontSize',18);
    
    % Create ylabel
    ylabel('Corresponding Dynamic Stiffness per unit surface [N/m^3]','FontSize',18);
    
    % Create xlabel
    xlabel('Force [N]','FontSize',18);
    
    %fig1 = create_mash(X3, Y3, Z3);
    %title(['Waterfall ',cell2mat([conf(i).conf.campione,...
    %    ' ',conf(i).conf.piastra,' ',conf(i).conf.punta])],'FontSize',18)
    % set(fig1.CurrentAxes ,'CLim',[7.14683052857818e-10 z_lim],'Colormap',...
    %     [0 0 0.515625;0 0 0.533622382198953;0 0 0.551619764397906;0 0 0.569617146596859;0 0 0.587614528795812;0 0 0.605611910994764;0 0 0.623609293193717;0 0 0.64160667539267;0 0 0.659604057591623;0 0 0.677601439790576;0 0 0.695598821989529;0 0 0.713596204188482;0 0 0.731593586387435;0 0 0.749590968586388;0 0 0.76758835078534;0 0 0.785585732984293;0 0 0.803583115183246;0 0 0.821580497382199;0 0 0.839577879581152;0 0 0.857575261780105;0 0 0.875572643979058;0 0 0.89357002617801;0 0 0.911567408376963;0 0 0.929564790575916;0 0 0.947562172774869;0 0 0.965559554973822;0 0 0.983556937172775;0 0.00155431937172756 1;0 0.0195517015706804 1;0 0.0375490837696333 1;0 0.0555464659685861 1;0 0.073543848167539 1;0 0.0915412303664919 1;0 0.109538612565445 1;0 0.127535994764398 1;0 0.14553337696335 1;0 0.163530759162303 1;0 0.181528141361256 1;0 0.199525523560209 1;0 0.217522905759162 1;0 0.235520287958115 1;0 0.253517670157068 1;0 0.27151505235602 1;0 0.289512434554973 1;0 0.307509816753926 1;0 0.325507198952879 1;0 0.343504581151832 1;0 0.361501963350785 1;0 0.379499345549738 1;0 0.39749672774869 1;0 0.415494109947643 1;0 0.433491492146596 1;0 0.451488874345549 1;0 0.469486256544502 1;0 0.487483638743455 1;0 0.505481020942408 1;0 0.523478403141361 1;0 0.541475785340314 1;0 0.559473167539267 1;0 0.57747054973822 1;0 0.595467931937173 1;0 0.613465314136125 1;0 0.631462696335078 1;0 0.649460078534031 1;0 0.667457460732984 1;0 0.685454842931937 1;0 0.70345222513089 1;0 0.721449607329843 1;0 0.739446989528796 1;0 0.757444371727749 1;0 0.775441753926702 1;0 0.793439136125655 1;0 0.811436518324608 1;0 0.829433900523561 1;0 0.847431282722514 1;0 0.865428664921467 1;0 0.88342604712042 1;0 0.901423429319373 1;0 0.919420811518326 1;0 0.937418193717279 1;0 0.955415575916232 1;0 0.973412958115185 1;0 0.991410340314138 1;0.00940772251309085 1 0.990592277486909;0.0274051047120438 1 0.972594895287956;0.0454024869109968 1 0.954597513089003;0.0633998691099498 1 0.93660013089005;0.0813972513089027 1 0.918602748691097;0.0993946335078557 1 0.900605366492144;0.117392015706809 1 0.882607984293191;0.135389397905762 1 0.864610602094238;0.153386780104715 1 0.846613219895285;0.171384162303668 1 0.828615837696332;0.189381544502621 1 0.810618455497379;0.207378926701574 1 0.792621073298426;0.225376308900527 1 0.774623691099473;0.243373691099479 1 0.756626308900521;0.261371073298432 1 0.738628926701568;0.279368455497385 1 0.720631544502615;0.297365837696338 1 0.702634162303662;0.315363219895291 1 0.684636780104709;0.333360602094244 1 0.666639397905756;0.351357984293197 1 0.648642015706803;0.36935536649215 1 0.63064463350785;0.387352748691103 1 0.612647251308897;0.405350130890056 1 0.594649869109944;0.423347513089009 1 0.576652486910991;0.441344895287962 1 0.558655104712038;0.459342277486915 1 0.540657722513085;0.477339659685868 1 0.522660340314132;0.495337041884821 1 0.504662958115179;0.513334424083774 1 0.486665575916226;0.531331806282727 1 0.468668193717273;0.549329188481679 1 0.450670811518321;0.567326570680632 1 0.432673429319368;0.585323952879585 1 0.414676047120415;0.603321335078538 1 0.396678664921462;0.62131871727749 1 0.37868128272251;0.639316099476443 1 0.360683900523557;0.657313481675396 1 0.342686518324604;0.675310863874349 1 0.324689136125651;0.693308246073301 1 0.306691753926699;0.711305628272254 1 0.288694371727746;0.729303010471207 1 0.270696989528793;0.74730039267016 1 0.25269960732984;0.765297774869112 1 0.234702225130888;0.783295157068065 1 0.216704842931935;0.801292539267018 1 0.198707460732982;0.819289921465971 1 0.180710078534029;0.837287303664923 1 0.162712696335077;0.855284685863876 1 0.144715314136124;0.873282068062829 1 0.126717931937171;0.891279450261782 1 0.108720549738218;0.909276832460734 1 0.0907231675392657;0.927274214659687 1 0.0727257853403129;0.94527159685864 1 0.0547284031413602;0.963268979057593 1 0.0367310209424074;0.981266361256545 1 0.0187336387434547;0.999263743455498 1 0.000736256544501934;1 0.982738874345549 0;1 0.964741492146596 0;1 0.946744109947644 0;1 0.928746727748691 0;1 0.910749345549738 0;1 0.892751963350785 0;1 0.874754581151833 0;1 0.85675719895288 0;1 0.838759816753927 0;1 0.820762434554974 0;1 0.802765052356022 0;1 0.784767670157069 0;1 0.766770287958116 0;1 0.748772905759163 0;1 0.730775523560211 0;1 0.712778141361258 0;1 0.694780759162305 0;1 0.676783376963352 0;1 0.6587859947644 0;1 0.640788612565447 0;1 0.622791230366494 0;1 0.604793848167541 0;1 0.586796465968589 0;1 0.568799083769636 0;1 0.550801701570683 0;1 0.53280431937173 0;1 0.514806937172778 0;1 0.496809554973825 0;1 0.478812172774872 0;1 0.460814790575919 0;1 0.442817408376967 0;1 0.424820026178014 0;1 0.406822643979061 0;1 0.388825261780108 0;1 0.370827879581156 0;1 0.352830497382203 0;1 0.33483311518325 0;1 0.316835732984297 0;1 0.298838350785345 0;1 0.280840968586392 0;1 0.262843586387439 0;1 0.244846204188486 0;1 0.226848821989534 0;1 0.208851439790581 0;1 0.190854057591628 0;1 0.172856675392675 0;1 0.154859293193723 0;1 0.13686191099477 0;1 0.118864528795817 0;1 0.100867146596864 0;1 0.0828697643979117 0;1 0.0625 0;1 0.0538461538461537 0;1 0.0451923076923073 0;1 0.036538461538461 0;1 0.0278846153846146 0;1 0.0192307692307683 0;1 0.0105769230769219 0;1 0.00192307692307558 0;0.993269230769229 0 0;0.984615384615383 0 0;0.975961538461537 0 0;0.96730769230769 0 0;0.958653846153844 0 0;0.949999999999998 0 0;0.941346153846151 0 0;0.932692307692305 0 0;0.924038461538458 0 0;0.915384615384612 0 0;0.906730769230766 0 0;0.898076923076919 0 0;0.889423076923073 0 0;0.880769230769227 0 0;0.87211538461538 0 0;0.863461538461534 0 0;0.854807692307688 0 0;0.846153846153841 0 0;0.837499999999995 0 0;0.828846153846149 0 0;0.820192307692302 0 0;0.811538461538456 0 0;0.80288461538461 0 0;0.794230769230763 0 0;0.785576923076917 0 0;0.776923076923071 0 0;0.768269230769224 0 0;0.759615384615378 0 0;0.750961538461532 0 0;0.742307692307685 0 0;0.733653846153839 0 0;0.724999999999993 0 0;0.716346153846146 0 0;0.7076923076923 0 0;0.699038461538454 0 0;0.690384615384607 0 0;0.681730769230761 0 0;0.673076923076914 0 0;0.664423076923068 0 0;0.655769230769222 0 0;0.647115384615375 0 0;0.638461538461529 0 0;0.629807692307683 0 0;0.621153846153836 0 0;0.61249999999999 0 0;0.603846153846144 0 0;0.595192307692297 0 0;0.586538461538451 0 0;0.577884615384605 0 0;0.569230769230758 0 0;0.560576923076912 0 0;0.551923076923066 0 0;0.543269230769219 0 0;0.534615384615373 0 0;0.525961538461527 0 0;0.51730769230768 0 0;0.508653846153834 0 0;0.5 0 0]);
    
    
    title(['Waterfall ',cell2mat([conf(i).conf.campione,...
        ' ',conf(i).conf.piastra,' ',conf(i).conf.punta])],'FontSize',18)
    ylim([0 y_max])
    zlim([0 z_lim])
    %xlim([0 5])
    set(gca, "CLim",([0 1]))
    set(gca,'Colormap',...
        [0 0 0.515625;0 0 0.533622382198953;0 0 0.551619764397906;0 0 0.569617146596859;0 0 0.587614528795812;0 0 0.605611910994764;0 0 0.623609293193717;0 0 0.64160667539267;0 0 0.659604057591623;0 0 0.677601439790576;0 0 0.695598821989529;0 0 0.713596204188482;0 0 0.731593586387435;0 0 0.749590968586388;0 0 0.76758835078534;0 0 0.785585732984293;0 0 0.803583115183246;0 0 0.821580497382199;0 0 0.839577879581152;0 0 0.857575261780105;0 0 0.875572643979058;0 0 0.89357002617801;0 0 0.911567408376963;0 0 0.929564790575916;0 0 0.947562172774869;0 0 0.965559554973822;0 0 0.983556937172775;0 0.00155431937172756 1;0 0.0195517015706804 1;0 0.0375490837696333 1;0 0.0555464659685861 1;0 0.073543848167539 1;0 0.0915412303664919 1;0 0.109538612565445 1;0 0.127535994764398 1;0 0.14553337696335 1;0 0.163530759162303 1;0 0.181528141361256 1;0 0.199525523560209 1;0 0.217522905759162 1;0 0.235520287958115 1;0 0.253517670157068 1;0 0.27151505235602 1;0 0.289512434554973 1;0 0.307509816753926 1;0 0.325507198952879 1;0 0.343504581151832 1;0 0.361501963350785 1;0 0.379499345549738 1;0 0.39749672774869 1;0 0.415494109947643 1;0 0.433491492146596 1;0 0.451488874345549 1;0 0.469486256544502 1;0 0.487483638743455 1;0 0.505481020942408 1;0 0.523478403141361 1;0 0.541475785340314 1;0 0.559473167539267 1;0 0.57747054973822 1;0 0.595467931937173 1;0 0.613465314136125 1;0 0.631462696335078 1;0 0.649460078534031 1;0 0.667457460732984 1;0 0.685454842931937 1;0 0.70345222513089 1;0 0.721449607329843 1;0 0.739446989528796 1;0 0.757444371727749 1;0 0.775441753926702 1;0 0.793439136125655 1;0 0.811436518324608 1;0 0.829433900523561 1;0 0.847431282722514 1;0 0.865428664921467 1;0 0.88342604712042 1;0 0.901423429319373 1;0 0.919420811518326 1;0 0.937418193717279 1;0 0.955415575916232 1;0 0.973412958115185 1;0 0.991410340314138 1;0.00940772251309085 1 0.990592277486909;0.0274051047120438 1 0.972594895287956;0.0454024869109968 1 0.954597513089003;0.0633998691099498 1 0.93660013089005;0.0813972513089027 1 0.918602748691097;0.0993946335078557 1 0.900605366492144;0.117392015706809 1 0.882607984293191;0.135389397905762 1 0.864610602094238;0.153386780104715 1 0.846613219895285;0.171384162303668 1 0.828615837696332;0.189381544502621 1 0.810618455497379;0.207378926701574 1 0.792621073298426;0.225376308900527 1 0.774623691099473;0.243373691099479 1 0.756626308900521;0.261371073298432 1 0.738628926701568;0.279368455497385 1 0.720631544502615;0.297365837696338 1 0.702634162303662;0.315363219895291 1 0.684636780104709;0.333360602094244 1 0.666639397905756;0.351357984293197 1 0.648642015706803;0.36935536649215 1 0.63064463350785;0.387352748691103 1 0.612647251308897;0.405350130890056 1 0.594649869109944;0.423347513089009 1 0.576652486910991;0.441344895287962 1 0.558655104712038;0.459342277486915 1 0.540657722513085;0.477339659685868 1 0.522660340314132;0.495337041884821 1 0.504662958115179;0.513334424083774 1 0.486665575916226;0.531331806282727 1 0.468668193717273;0.549329188481679 1 0.450670811518321;0.567326570680632 1 0.432673429319368;0.585323952879585 1 0.414676047120415;0.603321335078538 1 0.396678664921462;0.62131871727749 1 0.37868128272251;0.639316099476443 1 0.360683900523557;0.657313481675396 1 0.342686518324604;0.675310863874349 1 0.324689136125651;0.693308246073301 1 0.306691753926699;0.711305628272254 1 0.288694371727746;0.729303010471207 1 0.270696989528793;0.74730039267016 1 0.25269960732984;0.765297774869112 1 0.234702225130888;0.783295157068065 1 0.216704842931935;0.801292539267018 1 0.198707460732982;0.819289921465971 1 0.180710078534029;0.837287303664923 1 0.162712696335077;0.855284685863876 1 0.144715314136124;0.873282068062829 1 0.126717931937171;0.891279450261782 1 0.108720549738218;0.909276832460734 1 0.0907231675392657;0.927274214659687 1 0.0727257853403129;0.94527159685864 1 0.0547284031413602;0.963268979057593 1 0.0367310209424074;0.981266361256545 1 0.0187336387434547;0.999263743455498 1 0.000736256544501934;1 0.982738874345549 0;1 0.964741492146596 0;1 0.946744109947644 0;1 0.928746727748691 0;1 0.910749345549738 0;1 0.892751963350785 0;1 0.874754581151833 0;1 0.85675719895288 0;1 0.838759816753927 0;1 0.820762434554974 0;1 0.802765052356022 0;1 0.784767670157069 0;1 0.766770287958116 0;1 0.748772905759163 0;1 0.730775523560211 0;1 0.712778141361258 0;1 0.694780759162305 0;1 0.676783376963352 0;1 0.6587859947644 0;1 0.640788612565447 0;1 0.622791230366494 0;1 0.604793848167541 0;1 0.586796465968589 0;1 0.568799083769636 0;1 0.550801701570683 0;1 0.53280431937173 0;1 0.514806937172778 0;1 0.496809554973825 0;1 0.478812172774872 0;1 0.460814790575919 0;1 0.442817408376967 0;1 0.424820026178014 0;1 0.406822643979061 0;1 0.388825261780108 0;1 0.370827879581156 0;1 0.352830497382203 0;1 0.33483311518325 0;1 0.316835732984297 0;1 0.298838350785345 0;1 0.280840968586392 0;1 0.262843586387439 0;1 0.244846204188486 0;1 0.226848821989534 0;1 0.208851439790581 0;1 0.190854057591628 0;1 0.172856675392675 0;1 0.154859293193723 0;1 0.13686191099477 0;1 0.118864528795817 0;1 0.100867146596864 0;1 0.0828697643979117 0;1 0.0625 0;1 0.0538461538461537 0;1 0.0451923076923073 0;1 0.036538461538461 0;1 0.0278846153846146 0;1 0.0192307692307683 0;1 0.0105769230769219 0;1 0.00192307692307558 0;0.993269230769229 0 0;0.984615384615383 0 0;0.975961538461537 0 0;0.96730769230769 0 0;0.958653846153844 0 0;0.949999999999998 0 0;0.941346153846151 0 0;0.932692307692305 0 0;0.924038461538458 0 0;0.915384615384612 0 0;0.906730769230766 0 0;0.898076923076919 0 0;0.889423076923073 0 0;0.880769230769227 0 0;0.87211538461538 0 0;0.863461538461534 0 0;0.854807692307688 0 0;0.846153846153841 0 0;0.837499999999995 0 0;0.828846153846149 0 0;0.820192307692302 0 0;0.811538461538456 0 0;0.80288461538461 0 0;0.794230769230763 0 0;0.785576923076917 0 0;0.776923076923071 0 0;0.768269230769224 0 0;0.759615384615378 0 0;0.750961538461532 0 0;0.742307692307685 0 0;0.733653846153839 0 0;0.724999999999993 0 0;0.716346153846146 0 0;0.7076923076923 0 0;0.699038461538454 0 0;0.690384615384607 0 0;0.681730769230761 0 0;0.673076923076914 0 0;0.664423076923068 0 0;0.655769230769222 0 0;0.647115384615375 0 0;0.638461538461529 0 0;0.629807692307683 0 0;0.621153846153836 0 0;0.61249999999999 0 0;0.603846153846144 0 0;0.595192307692297 0 0;0.586538461538451 0 0;0.577884615384605 0 0;0.569230769230758 0 0;0.560576923076912 0 0;0.551923076923066 0 0;0.543269230769219 0 0;0.534615384615373 0 0;0.525961538461527 0 0;0.51730769230768 0 0;0.508653846153834 0 0;0.5 0 0]);
    colorbar(gca,'LineWidth',1);
    
    % impostazione visualizzazione isometrica
    set(gca,'View', [30 80])
    S.EdgeColor ='none';
    S.FaceColor ='interp';
    % salvataggio PDF figura isometrica
    exportgraphics(gcf,cell2mat(['Waterfall_SvsF-',conf(i).conf.punta,'_',...
        num2str(i),'_3d.pdf']),'BackgroundColor', 'none')
    
    % salvataggio formato fig
    savefig(cell2mat(['Waterfall_SvsF-',conf(i).conf.punta,'_',num2str(i),...
        '.fig'])); % cell2mat(['Collezione_Coerenza_',conf.campione,'_',conf.piastra]
    
    % impostazione visualizzazione dall'alto
    set(gca,'View', [0 90])
    S.EdgeColor ='none';
    S.FaceColor ='none';
    S.Marker = 'square';
    %S.MarkerSize
    S.MarkerFaceColor = 'Flat';
    % salvatoggio PDF figura dall'alto
    exportgraphics(gcf,cell2mat(['Waterfall_SvsF-',conf(i).conf.punta,'_',...
        num2str(i),'.pdf']), 'BackgroundColor', 'none')
    
    % salvataggio formato fig
    savefig(cell2mat(['Waterfall_SvsF-',conf(i).conf.punta,'_',num2str(i),...
        '.fig'])); % cell2mat(['Collezione_Coerenza_',conf.campione,'_',conf.piastra]
    
end
%%




%%
%<<<<<<<<<<<<<<
% Calibrazione
%<<<<<<<<<<<<<<
%massa di calibrazione
m_cal = 1.487; %#ok<NASGU>
%m_cal2 = piastre.massa('pesante2') + 0.0573;
m_cal = 2.951;
m_cal = piastre.massa( conf.conf.piastra ) ;

%Calcolo della massa dinamica normalizzata radice di PSD(F)/PSD(A) /
%massa di calibrazione

norm_mass = (PSD(1).F./PSD(1).A).^(1/2)/m_cal;
%norm_mass(:,4) = [];

figure, hold on
title (cell2mat(['Calibrazione - piastra ', conf(1).conf.piastra, ' - punta ', conf(1).conf.punta]), 'FontSize', 18)
xlabel ('Frequenza [Hz]', 'FontSize', 18);
ylabel ('Massa dinamica / massa sospesa', 'FontSize', 18);
%plot della massa dinamica normalizzata
plot (f, norm_mass,'Color', [0.5 0.5 0.5])
set(gca, 'XScale', 'log'), grid on
xlim ([10 1000]), ylim([0.5 1.2])

% Calcolo della massa dinamica media sulle ripetizioni
norm_mass_av = mean (norm_mass, 2);

% valore medio in frequenza della massa dinamica media
mediumval = mean (norm_mass_av(33:313));
%calib_m1 = 1./norm_mass_av

plot (f, norm_mass_av, 'k', 'LineWidth', 2)

plot (f, ones(size(f))*1.05, 'k-.', 'LineWidth', 1)
plot (f, ones(size(f)), 'k', 'LineWidth', 2)
plot (f, ones(size(f))*0.95, 'k-.', 'LineWidth', 1)

norm_A1 = PSD(1).A.*norm_mass_av.^2;
% norm_A2 = PSD(2).A.*norm_mass_av.^2;
% norm_mass2 = (PSD(2).F./norm_A2).^(1/2)/m_cal2;

%plot (f, norm_mass2, 'b', 'LineWidth', 1)
plot (f, (PSD(1).F./norm_A1).^(1/2)/m_cal, 'r', 'LineWidth', 1)
plot (f, norm_mass/mediumval,'b')
ylim ([0 1.2])
xlim ([80 1200])

calib = 1/mediumval;
save('calibrazione','calib')
%%
%<<<<<<<<<<<<<<<<<<<<<<<<
% Confronto in frequenza
%<<<<<<<<<<<<<<<<<<<<<<<<

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Confronto energia / massimo della forza
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure, hold on
for i=1:N
    [~,c] = size(sng(i).F_filt);
    energy = [];
    for j = 1:c
        temp = cumtrapz(sng(i).F_filt(:,j).^2 );
        energy(j) = temp(end);
    end
    
    
    plot(max(sng(i).F_filt), energy,'+', 'LineWidth', 3)
end
xlabel('Max force [N]')
ylabel('Energy')
title('Energy vs Force')
savefig('Energy_vs_Force.fig');


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
[PSD_F, ~]= periodogram(F, win_F, r, fs); %PSD Forza [N^2]
[PSD_F1, ~]= periodogram(F, win_1, r, fs); %PSD Forza [N^2]
[PSD_F2, ~]= periodogram(F_filt, win_F, r, fs); %PSD Forza [N^2]
[PSD_F3, f]= periodogram(F_filt, win_1, r, fs); %PSD Forza [N^2]

figure (56),hold on
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

%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% test della crosstrasformata per la media delle risposte
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

[pxy, ~] = cpsd(sng(1).F_filt, sng(1).A_filt, win_1, [], r, fs);

pxy_av = mean(pxy,2);

pxx = cpsd (sng(1).F_filt, sng(1).F_filt, win_1, [], r, fs);
pxx_av = mean (pxx,2);

k= (pxy_av).^(-1) .* ((1i*2*pi*f).^2);

txy = tfestimate(sng(1).F_filt,sng(1).A_filt, win_1, [], r, fs);
k_tfestimate = txy.^(1);

% figure (299)
% hold on
% plot (f, 20*log10(abs(k)))
% semilogx (f, 10 * log10(PSD(1).Kav),'color',string(colore(2,:)), 'LineWidth', 2)
%
% %plot della massa
% plot (f, 10*log10((piastre.massa(conf(i).conf.piastra)+0.0573)^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
% semilogx (f, 20 * log10(k_tfestimate.*(2*pi*f).^2),'color',string(colore(2,:)), 'LineWidth', 3)

figure (298)
hold on
plot (f, 20*log10 (abs (mean (txy, 2))), 'color', string(colore(3,:)), 'LineWidth', 2)
plot (f, 20*log10 (mean (abs (txy), 2)), '--', 'color', string(colore(3,:)), 'LineWidth', 2)
plot (f, 20*log10 (abs(pxy_av./pxx_av)), 'color', string(colore(1,:)), 'LineWidth', 2)
plot (f, 10*log10 (PSD(1).Aav./PSD(1).Fav), 'color', string(colore(2,:)), 'LineWidth', 2)
grid on
set(gca, 'XScale', 'log')
legend ('TFestimate','Cross transform k','PSD Kav')
xlabel (['Frequenza [Hz]'], 'FontSize', 18);
ylabel (['10 Log_{10} (K)'], 'FontSize', 18);

figure (299)
hold on
plot (f, -20*log10 (abs(mean (txy, 2)))+ 20*log10((2*pi*f).^2), 'color', string(colore(3,:)), 'LineWidth', 2)
plot (f, 20 *log10(abs(pxx_av./pxy_av) .* (2*pi*f).^2),'color', string(colore(1,:)), 'LineWidth', 2)
plot (f, 10 *log10(PSD(1).Kav),'color',string(colore(2,:)), 'LineWidth', 2)
%
%plot della massa
plot (f, 20*log10((piastre.massa(conf(i).conf.piastra))*(2*pi*f).^2),'b--', 'LineWidth', 2)
% semilogx (f, 20 * log10(k_tfestimate.*(2*pi*f).^2),'color',string(colore(2,:)), 'LineWidth', 3)

%plot della massa
%plot (f, 10*log10((piastre.massa(conf(i).conf.piastra)+0.0573)^2*(2*pi*f).^4),'b--', 'LineWidth', 2)
grid on
set(gca, 'XScale', 'log')
legend ('TFestimate','Cross transform k','PSD Kav')
xlabel (['Frequenza [Hz]'], 'FontSize', 18);
ylabel (['10 Log_{10} (K)'], 'FontSize', 18);

% Il risultato � molto simile tra tfestimate e cpsd. Per� cpsd sembra
% ottenere risultati migliori ad altissime frequenze dove alcune spike
% spurie vengono elise. La periodogram fornisce risultati simili a cpsd in
% termini di forma spettrale ma con un bias di circa 9dB. L'accelerazione
% in particolare, risulta pi� bassa se calcolata con cpsd e tfestimate.



%%

figure (301)
hold on
plot (f, 10*log10(abs(pxy(:,1))))
plot (f, 10*log10(PSD(1).A(:,1)))
set(gca, 'XScale', 'log')


%%
figure (300)
hold on
plot (f, 10*log10(abs(pxx(:,1) ./ pxy(:,1))))
plot (f, 10*log10(abs(pxx(:,1) ./ PSD(1).A(:,1))))
set(gca, 'XScale', 'log')
legend ('Cross transform','periodogram')


%%
% Calcolo ai sensi di quando proposto da:
% Pavement stiffness measurements in relation to mechanical impedance
set (0,'DefaultFigureWindowStyle','docked')

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
[campioni] = tabella_campioni (conf(1).conf,piastre);

N = length(result.indici_colpo.max_A);
fs = 52100;

[r,~] = size(sng(1).F(:,1));
dt = 1/fs; time1 = 1000*(0:dt:r/fs-dt);

%%integrazione accelerazione in velocit�
for i=1:N
    v_t = cumtrapz (sng(i).A)/fs;
    x = round(length(v_t)/2);
    y = v_t(x,:);
    v_drift = (ones(size(v_t)) .* y/x  );
    v_drift = cumtrapz (v_drift);
    
    %v_t = cumtrapz (sng(1).A_filt)/fs;
    v_t_corr = v_t - v_drift;
    
    
    
    colpo = 7;
    figure (700), hold on,
    yyaxis left
    plot (time1, sng(i).A_filt(:,colpo)),
    yyaxis right
    %plot(v_t(:,1))
    plot (time1, v_t_corr(:,colpo))
    xlim ([15 40])
    
    figure (701), hold on,
    plot (time1, sng(i).F_filt(:,colpo)),
    yyaxis right
    %plot(v_t(:,1))
    plot (time1, v_t_corr(:,colpo))
    xlim ([15 40])
    
    MI = max(sng(i).F) ./ max(v_t_corr);
    figure (702), hold on
    plot (10*log10(MI),'*')
    ylim ([ 0  1.2* max(10*log10( MI))])
    grid on
    
end
%%



v_0 = 0;
v_1 = v_0 + a_0 * 1/fs
v_2 = v1