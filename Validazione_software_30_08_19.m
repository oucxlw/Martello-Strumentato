%Martello strumentato 

set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables
load dati.mat

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
campione='segnale di test';
accelerometro=2;
punta='m';
piastra='p2';
martellamento='auto';

% x= p2_a2_m_mart;   % Force [N] 
% y= p2_a2_m_acc;    % Accelerazione [m/s^2]
% x=a2_a1_m_mart;
% y=a2_a1_m_acc2;

%<<<<<<<<<<<<<<<<<<<<<<<<
% Parametri di controllo
%<<<<<<<<<<<<<<<<<<<<<<<<

% Parametri fisici
sensit_F=1/1000;        % Sensibilità Martello strumentato [V/N]
sensit_A1=1/1000;%95.8/1000;    % Sensibilità Accelerometro1 [V/g]
sensit_A2=1/1000;%0.03/1000;   % Sensibilità Accelerometro2 [V/g]
g=9.81;                 % Accelerazione di gravità [m/s^2]
div_F=2000;             % Divisione applicata alla Forza prima della scrittura sul file WAVE
div_A=500;              % Divisione applicata alla Accelerazione prima della scrittura sul file WAVE

% Parametri di ricerca
fs=52100;               % Freq. di campionamento;
soglia=10;              % Soglia dei picchi;
delay=round(0.1*fs);    % Sample da saltare una volta superata la soglia
inizio=1;%15*fs;        % Punto di inizio della ricerca dei picchi;
fine=fs;%size(x);       % Punto di fine della ricerca dei picchi

% Parametri di filtro
bandwith0=0;            % Larghezza di banda richiesta al singolo colpo

% Dimensioni dei campioni
L_pre=round(1/2000*fs); % Lunghezza della parte prima del picco
L_coda=round(0.1*fs);   % Lunghezza della coda dei segnali

% Filtraggio doppi colpi
filt_doppi=1;           % Se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi

% Normalizzazione colpi
norm=0;                 % Se norm=1 i colpi vengono normalizzati
% Finestratura
finestratura=0;         % Tipo di finestratura da applicare
% 0 = nessuna finestratura
% 1 = finestratura quadrata
% 2 = hamming (da verificare)
% 3 = blackmanharris ritardata
% 4 = blackmanharris anticipata
% 5 = Hann Window

% Plotting
ascissamin=100;         % Frequenza minima da plottare nei grafici 
ascissamax=8000;        % Frequenza massima da plottare nei grafici
misura=['Campione ',num2str(campione),', ',martellamento,', punta ',punta,', piastra ',piastra,',Hann, PSDvsFFT, ',num2str(bandwith0),' Hz'];

%
% generazione segnale di prova
%


pos=1:fs;
force=normpdf(pos,300,35);
force=40*force/max(force);
figure (1)
subplot(2,1,1)
plot (pos,force)
L_force=length(force);
FFT_force=fft(force);
f = fs*(0:(L_force-1))/L_force;

k0=12e6;
FT=k0./((f*2*pi).^2);

FFT_acc = FT./FFT_force;
acc = ifft(1000*FFT_force);
figure (1)
hold on
subplot(2,1,1)
plot (pos,FFT_acc)
hold off

x=force';
y=acc';

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calibrazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% x=  div_F*x;     % Force [N] 
% y=g*div_A*(-y);

% L = length(x);
% dt=1/fs; time=[0:dt:L/fs-dt];
figure(100)
subplot(2,1,1), hold on, plot (x)
subplot(2,1,2), hold on, plot (y), 
hold off

%<<<<<<<<<<<<<<<<<<<<
% Ricerca dei PICCHI
%<<<<<<<<<<<<<<<<<<<<
picchi=[];
ii=inizio;
while   ii < fine(1)
    if  abs(x(ii)) > soglia
        picchi = [picchi; ii];
        ii=ii+delay;
    end
    ii=ii+1;
end
%Calcolo del numero di picchi
%picchi_s = picchi./fs;
n_picchi = length(picchi)

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
F=[]; A=[];
x_sel=[];y_sel=[]; % x_sel e y_sel contencono tutti i picchi selezionati in un unico vettore
scarti=0;
delay2=0;%round(2.3/1000*fs);
for j = 1: n_picchi
    in  = round(picchi(j) - L_pre);
    out = round(picchi(j) + L_coda);
    
    %<<<<<<<<<<<<<<<<<<<<<
    % Media Mobile
    %<<<<<<<<<<<<<<<<<<<<<
    x_temp = x(in:out);%movavg(x(in:out), 'exponential', 20);
    y_temp = y(delay2+in:delay2+out);    
  
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Filtraggio per doppio colpo
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    if filt_doppi==1
        n_colpi = findpeaks(x_temp,fs,'MinPeakDistance',0.002,'Threshold',0,'MinPeakHeight',5);
    else
        n_colpi = 1;  % I colpi non vengono filtrati
    end
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Filtraggio per intensità del colpo
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if  length(n_colpi) == 1 && (max(x_temp)>= 30) && (max(x_temp) <=70) % && (max(y_temp) >=50) && (max(y_temp) <=150) %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
        x_sel = [x_sel; x_temp];
        y_sel = [y_sel; y_temp];
        F = [F, x_temp];
        A = [A, y_temp];           
    else
%         figure(111) % Stampa dei colpi scartati
%         hold on
%         findpeaks(x_temp,fs,'MinPeakDistance',0.002,'Threshold',0,'MinPeakHeight',5);
%         hold off
scarti = scarti+1;
    end
%j = j+1;
end
picchi_sel = n_picchi - scarti

figure(110) % Stampa dei picchi selezionati
subplot(2,1,1), hold on, plot(x_sel)
subplot(2,1,2), hold on, plot(y_sel), 
hold off


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
F_filt=[]; A_filt=[]; V_filt=[];D_filt=[];      % Matrici Forza e Accelerazione picchi selezionati
F_filt_all=[];  A_filt_all=[]; V_filt_all=[]; D_filt_all=[];  % Vettori con picchi selezionati in serie (serve per calcolare la coerenza)
F_norm=[];  A_norm=[]; V_norm=[]; D_norm=[];      % Matrici Forza e Accelerazione normalizzate

L = length(F(:,1));
dt=1/fs; time=1000*(0:dt:L/fs-dt);

for kk=1:picchi_sel
    cont1=0;  
  
    %sottrazione della media
    F(:,kk)=F(:,kk)-mean(F(:,kk));
    A(:,kk)=A(:,kk)-mean(A(:,kk));
    
    %<<<<<<<<<<<<<<<<<<<
    % Finestratura di F
    %<<<<<<<<<<<<<<<<<<<
    switch finestratura
        case 0
            F_filt = [F_filt,F(:, kk)];
        
        case 1 %FILTRO: ad onda quadra
            Lw1=length(F(:,kk));
            [m_x, m_y]=find(F(:,kk) == max(F(:,kk)));
            for i=1:m_x
                if F(i,kk)<0.001
                    cont1=cont1+1;
                end
            end
            z=0;
            while z==0        
                for i=m_x:Lw1        
                    if  F(i,kk)<0
                        z=i; break;
                    end
                end
            end
            %w1=[]; 
            Liniz=cont1; Lwind=z-1-cont1;
            w1=[zeros(Liniz, 1); ones(Lwind, 1);  zeros(Lw1-(Lwind+Liniz), 1)];  
            F_filt = (F(:, kk)) .* w1;    
            
        case 2 % finestra hamming
            Lw1=length(F(:,kk));
            [m_x, m_y]=find(F(:,kk) == max(F(:,kk)));
            for i=1:m_x
                if F(i,kk)<0.001
                    cont1=cont1+1;
                end
            end
            z=0;
            while z==0        
                for i=m_x:Lw1        
                    if  F(i,kk)<0
                        z=i; break;
                    end
                end
            end
            Liniz=cont1; Lwind=z-1-cont1;
            H=hamming(Lwind);         
            figure (200)
            plot(H)
            for ii=1:Lw1
                 if ii>Liniz && ii< (Liniz+Lwind)
                     F_filt=[F_filt,F(ii:kk).*H(ii-Liniz)];
                 else
                     F_filt(ii)=0;
                 end
            end
            figure (202)
            plot(F_filt)
            
        case 3 %bl_harr caso 1 
            t0=L_pre-3; % inizio della finestra di salita
            [~,pos] = max(F(:,kk));
            jj=pos+1;
            while F(jj,kk)>0
                jj=jj+1;
            end
            tf = jj;
            Lp = tf - t0;
            Delta_f = round(Lp/10);
            Sgn=F(1:tf+Delta_f-1,kk);
            win0= blackmanharris(2*t0);
            win0=win0(1:t0);
            uno = ones(Lp-2,1);
            winf= blackmanharris(Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn=Sgn.*win;
            zeri=zeros(length(F(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
                F_filt=[F_filt, temp];
            
        case 4 %bl_harr caso 2 
            t0=L_pre-3; % inizio della finestra di salita
            [~,pos] = max(F(:,kk));
            jj=pos+1;
            while F(jj,kk)>0
                jj=jj+1;
            end
            tf = jj;
            Lp = tf - t0;
            Delta_f = round(Lp/5);
            Sgn=F(1:tf-1,kk);
            win0= blackmanharris(2*t0);
            win0=win0(1:t0);
            uno = ones(Lp-2-Delta_f,1);
            winf= blackmanharris(Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn=Sgn.*win;
            zeri=zeros(length(F(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            F_filt=[F_filt, temp];
                
        case 5 %hann
            t0=L_pre-3; % inizio della finestra di salita
            [~,pos] = max(F(:,kk));
            jj=pos+1;
            while F(jj,kk)>0
                jj=jj+1;
            end
            tf = jj;
            Lp = tf - t0;
            Delta_f = round(Lp/5);
            Sgn=F(1:tf-1,kk);
            win0= hann(2*t0);
            win0=win0(1:t0);
            uno = ones(Lp-2-Delta_f,1);
            winf= hann (Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn=Sgn.*win;
            zeri=zeros(length(F(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            F_filt=[F_filt, temp];

        otherwise
          F_filt =[F_filt, F(:, kk)];   
    end
%    figure
%    plot (F_filt(:,kk)), hold on, plot(F(:,kk),'k'), xlim([0 500])   
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Calcolo NORMA F per il singolo colpo
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if norm==1
        temp=0;
        for ll=1:length(F(:,kk))
            temp=temp+(F(ll,kk)^2);        
        end
        norma_F = temp^(1/2);
    else
        norma_F = 1;
    end
    F_filt(:,kk) = F_filt(:,kk)./norma_F(end);
    F_filt_all = [F_filt_all; F_filt(:,kk)]; 

    %<<<<<<<<<<<<<<<<<<<
    % Finestratura di A
    %<<<<<<<<<<<<<<<<<<<
    switch finestratura
        case 0 % Exponenziale:           
            A_filt= [A_filt, A(:, kk)];      
            A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A(:,kk)];
            
            V_filt=[V_filt, cumtrapz(time,A(:,kk))];       %Calcolo Velocità V a partire da A non finestrato
            D_filt=[D_filt, cumtrapz(time,V_filt(:,kk))];    %Calcolo Spostamento D a partire da V non finestrato
            V_filt(:,kk) = V_filt(:,kk)/norma_F(end);       %NON finestro e nel caso normalizzo V
            D_filt(:,kk) = D_filt(:,kk)/norma_F(end);       %NON finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];       
       
        case 4 %bl_harr caso 2 
            t0=L_pre-3; % inizio della finestra di salita
            [~,pos] = max(A(:,kk));
            jj=pos+1;
            while A(jj,kk)>0
                jj=jj+1;
            end
            tf = jj;
            Lp = tf - t0;
            Delta_f = round(Lp/2);
            Sgn=A(1:tf-1,kk);
            win0= blackmanharris(2*t0);
            win0=win0(1:t0);
            uno = ones(Lp-2-Delta_f,1);
            winf= blackmanharris(Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn=Sgn.*win;
            zeri=zeros(length(A(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            A_filt=[A_filt, temp];     
            A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A_filt(:,kk)];
            
            V_filt=[V_filt, cumtrap(time,A(:,kk))];       %Calcolo Velocità V a partire da A non finestrato
            D_filt=[D_filt, cumtrap((1/fs),V_filt(:,kk))];    %Calcolo Spostamento D a partire da V non finestrato
            V_filt(:,kk) = win.*V_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo V
            D_filt(:,kk) = win.*D_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];          %Accodo D
            
        case 5 %hann accelerazione 
            t0=L_pre-3; % inizio della finestra di salita
            [~,pos] = max(A(:,kk));
            jj=pos+1;
            while A(jj,kk)>0
                jj=jj+1;
            end
            tf = jj;
            Lp = tf - t0;
            Delta_f = round(Lp/2);
            Sgn_A=A(:,kk);
            win0= hann(2*t0);
            win0=win0(1:t0);
            uno = ones(length(Sgn_A)-t0-Delta_f-1,1);
            winf= hann(Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn_A=Sgn_A.*win;
            A_filt=[A_filt, Sgn_A]; %Accelerazione
            A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A_filt(:,kk)];      
            Lung=length(A_filt(:,1));
            f =0:fs/Lung:fs/2;

            V_filt=[V_filt, cumtrapz(time,A(:,kk))];        %Calcolo Velocità V a partire da A non finestrato
            D_filt=[D_filt, cumtrapz(time,V_filt(:,kk))];   %Calcolo Spostamento D a partire da V non finestrato
            V_filt(:,kk) = win.*V_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo V
            D_filt(:,kk) = win.*D_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];          %Accodo D
            
        otherwise %NESSUNA FINESTRATURA
            A_filt= [A_filt, A(:, kk)];      
            A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A(:,kk)];
            
            V_filt=[V_filt, cumtrapz(time,A(:,kk))];       %Calcolo Velocità V a partire da A non finestrato
            D_filt=[D_filt, cumtrapz(time,V_filt(kk))];    %Calcolo Spostamento D a partire da V non finestrato
            V_filt(:,kk) = V_filt(:,kk)/norma_F(end);       %NON finestro e nel caso normalizzo V
            D_filt(:,kk) = D_filt(:,kk)/norma_F(end);       %NON finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];          %Accodo D
            
    end
    
    if kk == 1 && finestratura~=0
        figure (1),hold on
        subplot(2,1,1), plot (F(:,kk)), hold on, plot(F_filt(:,kk), 'r-.', 'LineWidth', 2), plot(win.*40), xlim([0 400]), grid on,
        subplot(2,1,2), plot(A(:,kk)), hold on, plot(A_filt(:,kk), 'r-.', 'LineWidth', 2), plot(win.*100),xlim([0 1000]), grid on,
    end
    %kk=kk+1;
end
hold off
%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle FRFs
%<<<<<<<<<<<<<<<<<<<<
PSD_F=[];   PSD_A=[];   PSD_V=[];   PSD_D=[]; 
FFT_F=[];   FFT_A=[];   FFT_V=[];   FFT_D=[];

Ffilt_all2=[];
Afilt_all2=[];

L = length(F_filt(:,1));
dt=1/fs; time=1000*(0:dt:L/fs-dt);

scarti2=0;
for j=1:picchi_sel 

    %Spettro di densità di potenza (PSD vs. freq) di forza e accelerazione:
    [Pxx, f] = periodogram(F_filt(:,j), [], L, fs); %PSD Forza [N^2]
           
    f0=find(f>ascissamin,1);
    fmax=find(Pxx(f0:end)<((Pxx(f0)/10)),1); 
    fmax=f(fmax+f0);
    
    %??????????????????????????????????????????????????????????????????????????????????????????
    % Sto cercando la larghezza di banda con la condizione:`find(Pxx(f0:end)<(Pxx(f0)/10),1)`
    % Pxx però è un psd quindi quadratica. Non sarebbe meglio utilizzare:
    % find(sqrt(Pxx(f0:end))<(sqrt(Pxx(f0))/10),1)
    %??????????????????????????????????????????????????????????????????????????????????????????
    
    if  fmax>bandwith0   %se la richiesta di banda è soddisfatta calcolo gli spettri e li colleziono

        %Calcolo e memorizzo Spettri in fft di F, A,V e D
        FFT_F=[FFT_F, fft(F_filt(:,j))];
        FFT_A=[FFT_A, fft(A_filt(:,j))]; %FFT accelerazione
        FFT_V=[FFT_V, fft(V_filt(:,j))]; %FFT velocità 
        FFT_D=[FFT_D, fft(D_filt(:,j))]; %FFT spostamento=Displacement
        
        %Calcolo Spettri di densità di potenza (PSD vs. freq)
        % Periodogram restituisce il modulo quadro della PSD del segnale
        % per unità di frequenza.
        [Pyy1, f] = periodogram(A_filt(:,j), [], L, fs); %PSD Accelerazione[(m/s^2)^2=m^2/s^4]
        [Pyy2, f] = periodogram(V_filt(:,j), [], L, fs); %PSD Velocità [(m/s)^2=m^2/s^2]
        [Pyy3, f] = periodogram(D_filt(:,j), [], L, fs); %PSD Spostamento =Displacement [m^2]
        
        %Memorizzo i PSD di F, A, V e D
        PSD_F = [PSD_F, Pxx]; 
        PSD_A = [PSD_A, Pyy1]; %PSD Accelerazione
        PSD_V = [PSD_V, Pyy2]; %PSD Velocità
        PSD_D = [PSD_D, Pyy3]; %PSD Spostamento=Displacement
        
        figure(101), grid on,
        sgtitle(misura)
        subplot(2,2,1), hold on, plot(time, F_filt(:, j-scarti2)), xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Force'), grid on, xlim([0 10])
        subplot(2,2,3), hold on, semilogx (f, 10*log10(PSD_F(:,j-scarti2))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 N/Hz)'), title('PSD Force'), 
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        subplot(2,2,2), hold on, plot(time, A_filt(:, j-scarti2)), xlabel('Time [ms]'), ylabel('Amplitude [g]'), title('Acceleration'), grid on, xlim([0 10])
        subplot(2,2,4), hold on, semilogx (f, 10*log10(PSD_A(:,j-scarti2))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration'), 
        grid on, set(gca, 'XScale', 'log'), xlim([10 20000])

        Ffilt_all2=[Ffilt_all2;F_filt(:,j)];
        Afilt_all2=[Afilt_all2;A_filt(:,j)];
        
    else        %Altrimeni aumento il contatore scarti2
        scarti2=scarti2+1;
    end
    %j=j+1;
end
hold off
%tolgo dal conteggio gli scarti in base alla larghezza di banda
picchi_sel = picchi_sel - scarti2

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% calcolo la media degli spettri
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%PSD
PSD_Fav = mean(sqrt(PSD_F),2)';
PSD_Aav = mean(sqrt(PSD_A),2)';
PSD_Vav = mean(sqrt(PSD_V),2)';
PSD_Dav = mean(sqrt(PSD_D),2)';
%FFT
FFT_Fav=mean(FFT_F,2)'; %forza
FFT_Aav=mean(FFT_A,2)'; %accelerazione
FFT_Vav=mean(FFT_V,2)'; %forza
FFT_Dav=mean(FFT_D,2)'; %accelerazione

%Plot degli spettri medi in PSD
figure (101), subplot(2,2,3), hold on, plot (f, 20*log10(PSD_Fav), 'b-.', 'LineWidth', 3)
figure (101), subplot(2,2,4), hold on, plot (f, 20*log10(PSD_Aav), 'b', 'LineWidth', 3)
saveas (gcf, ['Segnali e spettri-C_',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.fig'])

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo della frequenza massima
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
PSD_Fav_dB = 20.*log10(PSD_Fav);
fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));

fmax=f(fmax(1)+f0);
%plot sulla PSD della forza

figure (101), subplot (2,2,3),hold on

xl=xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'])
xl.LabelVerticalAlignment = 'bottom';
hold off

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% COERENZA usando Forza / Accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[Cxy,f] = mscohere(Ffilt_all2, Afilt_all2, round(length(Ffilt_all2)./picchi_sel),[],L,fs);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calcolo DYNAMIC MASS Mechanical Impedance Dynamic Stiffness
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Dynamic Mass
DMASS_av = PSD_Fav./PSD_Aav; %Modulo della Dynamic Mass;
save (['DMASS_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.mat'], 'DMASS_av');
DMASS_av_ph = angle(FFT_Fav./FFT_Aav); %trovo la fase media usando le FFT;

% Mechanical Impedance
MI_av = PSD_Fav./PSD_Vav; %Modulo dell'Impadenza meccanica
save (['MI_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.mat'], 'MI_av');
MI_av_ph = angle(FFT_Fav./FFT_Vav); %trovo la fase media usando le FFT;


% Dynamic Stiffness
Dstiff_av = PSD_Fav./PSD_Dav; %Modulo della Dynamic Stiffness
save (['Dstiffness_av, misura-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.mat'], 'Dstiff_av');
Dstiff_av_ph = angle(FFT_Fav./FFT_Dav); %trovo la fase media usando le FFT;


%<<<<<<<<<<<<<<<<<<<
% Plot Dynamic Mass
%<<<<<<<<<<<<<<<<<<<

clf
figure(104), 
sgtitle(misura)
subplot(3,1,1), plot(f,Cxy), set(gca, 'XScale', 'log'), 
title('Magnitude-Squared Coherence')
xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 

figure (104), subplot(3,1,2), hold on, semilogx (f, 20*log10(DMASS_av), 'LineWidth', 3),

k=12.5e6;

semilogx (f, 20*log10(k./((2*pi*f').*(2*pi*f'))), 'LineWidth', 3),
set(gca, 'XScale', 'log'), 
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Mass| (dB ref 1 kg)'), title(['Dynamic Mass (Force/Acceleration) Amplitude']), 
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']), grid on, xlim([ascissamin ascissamax])

figure (104), subplot(3,1,3), 
hold on,
plot (f, 180.*DMASS_av_ph(1:L/2+1)./pi, 'LineWidth', 3),
set(gca, 'XScale', 'log')
xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Mass (Force/Acceleration) Phase']),
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'])
grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
saveas (gcf, ['Coerenza e dmass-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.fig'])

E=k*(0.045)/(pi*0.01^2)/10^6


%<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Mechanical Impedence
%<<<<<<<<<<<<<<<<<<<<<<<<<<<
figure(105), 
sgtitle(misura)
subplot(3,1,1), plot(f,Cxy), set(gca, 'XScale', 'log'), 
title('Magnitude-Squared Coherence')
xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 

figure (105), subplot(3,1,2), hold on, semilogx (f, 20*log10(DMASS_av.*(2*pi*f')), 'LineWidth', 3), 
%semilogx (f, 20*log10(MI_av), 'LineWidth', 3),
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylim([0 100])
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Mech. Impedance| (dB ref 1 N s/m]'), title(['Mechanical Impedance (Force/Velocity) Amplitude']), 
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']), grid on, xlim([ascissamin ascissamax])

figure (105), subplot(3,1,3), 
hold on, plot (f, 180.*MI_av_ph(1:L/2+1)./pi, 'LineWidth', 3),
set(gca, 'XScale', 'log')
xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Mechanical Impedance (Force/Velocity) Phase']),
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'])
grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
saveas (gcf, ['Coerenza e MI-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.fig'])

%<<<<<<<<<<<<<<<<<<<<<<<<
% Plot Dynamic Stiffness
%<<<<<<<<<<<<<<<<<<<<<<<<
figure(106), 
sgtitle(misura)
subplot(3,1,1), plot(f,Cxy), set(gca, 'XScale', 'log'), 
title('Magnitude-Squared Coherence')
xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 

figure (106), subplot(3,1,2), hold on, semilogx (f, 20*log10(DMASS_av.*(2*pi*f').^2), 'LineWidth', 3),

%semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylim([100 150])
xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'), title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']), grid on, xlim([ascissamin ascissamax])

figure (106), subplot(3,1,3), 
hold on,
plot (f, 180.*Dstiff_av_ph(1:L/2+1)./pi, 'LineWidth', 3),
set(gca, 'XScale', 'log')
xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Stiffness (Force/Displacement) Phase']),
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'])
grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
saveas (gcf, ['Coerenza e Dstiff-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.fig'])                    


picchi_sel

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
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']), grid on, xlim([ascissamin ascissamax])
legend('Media sulle PSD')
figure (106), subplot(2,1,2), 
hold on,
plot (f, 180.*Dstiff_av_ph(1:L/2+1)./pi, 'LineWidth', 3),
set(gca, 'XScale', 'log')
xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Stiffness (Force/Displacement) Phase']),
xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'])
grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
%saveas (gcf, ['Coerenza e Dstiff-C',num2str(campione),'-',martellamento,'-',punta,'-',piastra,'-blackmanharris 2-PSDvsFFT-',num2str(bandwith0),'Hz','.fig'])                    






