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

[RR,CC]=size(F);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%Faccio un calcolo di F_filt per ottenere L_win 
[F_filt, L_win] = finestra_forza (F, window_F, fs);

%Finestro sia Accelerazione che Forza utilizzando finestra_accelerazione 5
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
picchi_sel2 = picchi_sel1 - scarti
PSD_A(:,tagli)=[];

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Filtraggio in intensità PSD
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%calcolo del massimo e del minimo dei plateaux delle PSD in N (quindi
%operando la radice)
bin=round(sqrt(picchi_sel2))+1;
Max_pic=sqrt(max(max(PSD_F)));
Min_pic=sqrt(min(max(PSD_F)));
delta=(Max_pic-Min_pic)/bin;
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

[Y,E] = discretize(sqrt(max(PSD_F)),bin);
values=1:bin;

figure,histfit(sqrt(max(PSD_F)),bin);


kkk=0;indice=0;
[r,c]=size(PSD_F);

for indice = 1:bin
    kkk=kkk+1;
    PSD_Fbin=[];
    PSD_Abin=[];
    for jj=1:c
        if Y(jj)==indice
            PSD_Fbin=[PSD_Fbin,PSD_F(:,jj)];
            PSD_Abin=[PSD_Abin,PSD_A(:,jj)];
        end
    end
  
    [RR,CC]=size(PSD_Fbin);
    
%     % Stampo dei picchi selezionati
%     x_sel = reshape(F,[],1);
%     y_sel = reshape(A,[],1);
%     l = length(x_sel);
%     dt=1/fs; time=1000*(0:dt:l/fs-dt);
%     figure(2) 
%     subplot(2,1,1), hold on, plot(time, x_sel)
%     subplot(2,1,2), hold on, plot(time, y_sel), 
%     hold off

    if CC>=1 & E(indice)>0.01
%         % Matrici F e An filtrate con il bandwidth
%         F_filt2 = F_filt ;
%         A_filt2 = A_filt;
%         F_filt2(:,tagli)=[];
%         A_filt2(:,tagli)=[];

        %Calcolo e memorizzo Spettri in fft di F, A, V e D
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % calcolo la media degli spettri
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %PSD
        PSD_Fav = mean(sqrt(PSD_Fbin), 2);
        PSD_Aav = mean(sqrt(PSD_Abin), 2);
        PSD_V1av = PSD_Aav./(2*pi*f); %velocità
        PSD_D1av = PSD_V1av./(2*pi*f); % displacement
        
%         %FFT
%         FFT_Fav = mean( fft (F_filt2,  [], 1), 2); %FFT forza
%         FFT_Aav = mean( fft (A_filt2, [], 1), 2); %FFT accelerazione
%         FFT_V1av = FFT_Aav(1:L/2+1)./(1i*2*pi*f); %velocità
%         FFT_D1av=FFT_V1av(1:L/2+1)./(1i*2*pi*f); % displacement

        %Ciclo for per plottare segnali e spettri (PSD):
        dt=1/fs; time1=1000*(0:dt:L/fs-dt);

        figure(101), grid on,
        sgtitle(misura)
        for j=1:c
            if Y(j)==indice
                subplot(2,2,1), hold on,
                plot(time1, F_filt(:, j),'color',string(colore(kkk,:))), xlabel('Time [ms]'), ylabel('Amplitude [N]'), title('Force'),
                grid on, xlim([0 10])
                
                subplot(2,2,3), hold on,
                semilogx (f, 10*log10(PSD_F(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 N/Hz)'), title('PSD Force'), 
                grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
                
                subplot(2,2,2), hold on,
                plot(time1, A_filt(:, j),'color',string(colore(kkk,:))), xlabel('Time [ms]'), ylabel('Amplitude [g]'), title('Acceleration'),
                grid on, xlim([0 10])
                
                subplot(2,2,4), hold on,
                semilogx (f, 10*log10(PSD_A(:, j))), xlabel('log(Frequency) [Hz]'), ylabel('20 log |PSD| (dB ref 1 m/s^2 Hz)'), title('PSD Acceleration'), 
                grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
            end
        end
        hold off

        %Plot degli spettri medi in PSD
        figure (101), subplot(2,2,3), hold on,
        plot (f, 20*log10(PSD_Fav), '-.','color',string(colore(kkk,:)), 'LineWidth', 3),grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        figure (101), subplot(2,2,4), hold on,
        plot (f, 20*log10(PSD_Aav), 'b' ,'color',string(colore(kkk,:)), 'LineWidth', 3),grid on, set(gca, 'XScale', 'log'), xlim([10 20000])
        
        saveas (gcf, ['Segnali e spettri-C_',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.fig'])

        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Calcolo della frequenza massima
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        PSD_Fav_dB = 20.*log10(PSD_Fav);
        fmax=find(PSD_Fav_dB(f0:end) <(PSD_Fav_dB(f0)-10));
        fmax=f(fmax(1)+f0);
        %plot sulla PSD della forza
        figure (101), subplot (2,2,3),hold on
        xl=xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
        hold off
        
%         %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%         % COERENZA usando Forza / Accelerazione
%         %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%         F_filtall = reshape(F_filt2, [],1);
%         A_filtall = reshape(A_filt2, [],1);
%         [r,c]=size(F_filt2);
%         [Cxy1,f] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],L,fs);
%         save ([num2str(round(indice)),'Coherence, misura-C',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Cxy1');
% 
%         clear F_filt2
%         clear A_filt2

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


        %_____________________________________________________
        % Dynamic Stiffness
        Dstiff1_av = PSD_Fav./PSD_D1av; %Modulo della Dynamic Stiffness
        save (['Dstiffness_av, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av');
%         Dstiff1_av_ph = angle(FFT_Fav(1:L/2+1)./FFT_D1av); %trovo la fase media usando le FFT;
%         save (['Dstiffness_av_ph, misura - ',num2str(campione),'-Acc',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'-',num2str(bandwidth),'Hz','.mat'], 'Dstiff1_av_ph');

        %<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot Dynamic Stiffness
        %<<<<<<<<<<<<<<<<<<<<<<<<
%         figure(106), 
%         sgtitle(misura)
%         subplot(3,1,1),
%         plot(f,Cxy1),
%         set(gca, 'XScale', 'log'), title('Magnitude-Squared Coherence')
%         xlabel('log(Frequency) [Hz]'), ylabel('[-]'), 
%         grid on, ylim([0 1.1]), xlim([ascissamin ascissamax]) 
% 
%         figure (106), %subplot(3,1,2),
%         hold on,
%         semilogx (f, 20*log10(Dstiff1_av),'color',string(colore(kkk,:)), 'LineWidth', 3),
%         %semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
%         set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylim([100 200])
%         xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'), title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
%         xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']); grid on, xlim([ascissamin ascissamax])

%         figure (106), subplot(3,1,3), 
%         hold on,
%         plot (f, 180.*Dstiff1_av_ph./pi,'color',string(colore(kkk,:)), 'LineWidth', 3),
%         set(gca, 'XScale', 'log')
%         xlabel('log(Frequency) [Hz]'), ylabel('Phase [°]'), title(['Dynamic Stiffness (Force/Displacement) Phase']),
%         xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz']);
%         grid on, xlim([ascissamin ascissamax]), ylim([-180 180]), hold off
%         saveas (gcf, ['Coerenza e Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,', ',num2str(bandwidth),'Hz','.fig'])                    

        %plot singole ripetizioni nel range selezionato
        
        figure (round(indice+200)),hold on
        for iii=1:CC
        semilogx (f, 10*log10( PSD_Fbin(:,iii).*((2*pi*f).^4)./PSD_A (:,iii) ),'color',string(colore(kkk,:)), 'LineWidth', 1),
        %semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
        end
        set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
        titolo=['Dynamic Stiffness (Force/Displacement) Amplitude (fron ',num2str(E(indice)),' to ',num2str(E(indice+1)),' N)'];
        title([titolo]),
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        grid on, xlim([20 ascissamax]),ylim([120 190])
        saveas (gcf, ['fron ',num2str(E(indice)),' to ',num2str(E(indice+1)),' N)',' N Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'.fig'])
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<
        % Plot singolo D-Stiffness
        %<<<<<<<<<<<<<<<<<<<<<<<<<<

        figure (107),hold on,
        semilogx (f, 20*log10(Dstiff1_av),'color',string(colore(kkk,:)), 'LineWidth', 3),
        %semilogx (f, 20*log10(Dstiff_av), 'LineWidth', 3), 
        set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
        xlabel('log(Frequency) [Hz]'), ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'), 
        title(['Dynamic Stiffness (Force/Displacement) Amplitude, Sample: ',campione,'']), 
        xline(fmax,'.',['Limite in frequenza: ',num2str(round(fmax)),' Hz'],'color',string(colore(kkk,:)));
        

    end
    
    figure (107)     
    grid on, xlim([20 ascissamax]),ylim([120 190])
    saveas (gcf, ['Collezione Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'.fig'])
    saveas (gcf, ['Collezione Dstiff-C',num2str(campione),'-Acc_',num2str(accelerometro),'-',martellamento,'-',punta,'-',piastra,'.png'])

end
