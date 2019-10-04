% Martello strumentato filtro intensita
set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables
load dati.mat

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
campione='teflon senza biadesivo, rip 1';
accelerometro=0;
punta='m'; %m=metallica; p=plastica; g=gomma.
piastra='pesante';
martellamento='auto';

x = pp_m_teflon_1 (:,1); % Force [N] 
y = pp_m_teflon_1 (:,accelerometro+2); % Accelerazione [m/s^2] 

%x = polipropilene_piccola_m_biad_3 (:,1); % Force [N] 
%y = polipropilene_piccola_m_biad_3 (:,accelerometro+2); % Accelerazione [m/s^2] 

%  x = slab_piccola_m_res_4 (:,1); % Force [N] 
%  y = slab_piccola_m_res_4 (:,accelerometro+2); % Accelerazione [m/s^2] 
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
inizio=1*fs;          % Punto di inizio della ricerca dei picchi;
fine=round(0.95*size(x));           % Punto di fine della ricerca dei picchi
% Parametri di filtro
bandwidth=0;            % Larghezza di banda richiesta al singolo colpo
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
ascissamin=10;         % Frequenza minima da plottare nei grafici 
ascissamax=4000;        % Frequenza massima da plottare nei grafici
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
