function [F_filt, L_win] = finestra_forza (Forze, window_F, fs)
%Function che esegue la finestrature delle colonne delle matrici Forza e Accelerazioni,
%permettendo di scegliere diverse tra varie finestre.
%Prende in Input: 
%I1) la matrice delle Forze (F);
%I2) il tipo di finestra selezionata per le Forze (window_F=numero tra 1 e 5), scegliendo tra: 
        %case 1 = onda quadra come finestra,
        %case 2 = funzione di hamming come finestra,
        %case 3 = funzione di Blackman-Harris (versione 1) come finestra,
        %case 4 = funzione di Blackman-Harris (versione 2) come finestra,
        %case 5 = funzione di hann come finestra;
%I3) La frequenza di campionamento del segnale in Hz
        
%Restituendo come Output: 
%O1) la matrice Forze dopo la finestratura (F_filt);
%O2) la lunghezza del filtro usato;


F=Forze;

%<<<<<<<<<<<<<<
% Finestratura
%<<<<<<<<<<<<<<
F_filt=[];       % Matrice Forza dei picchi selezionati

L = length(F(:,1));
dt=1/fs; time=1000*(0:dt:L/fs-dt);
[r,picchi_sel]=size(F);
for kk=1:picchi_sel
    cont1=0;  
  
    %sottrazione della media
    F(:,kk)=F(:,kk)-mean(F(:,kk));
    
    %<<<<<<<<<<<<<<<<<<<
    % Finestratura di F
    %<<<<<<<<<<<<<<<<<<<
    switch window_F    
        case 1 % finestra: onda quadra
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
            
        case 2 % finestra: hamming
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
%             figure (200)
%             plot(H)
            for ii=1:Lw1
                 if ii>Liniz && ii< (Liniz+Lwind)
                     F_filt=[F_filt,F(ii:kk).*H(ii-Liniz)];
                 else
                     F_filt(ii)=0;
                 end
            end
%             figure (202)
%             plot(F_filt)
            
        case 3 % finestra: Blackman-Harris (versione 1)
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
            
        case 4  %finestra: Blackman-Harris (versione 2)
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
                
        case 5 % finestra: hann
            t0=1;%L_pre-3; % inizio della finestra di salita
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
            L_win=length(win);
            Sgn=Sgn.*win;
            zeri=zeros(length(F(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            F_filt=[F_filt, temp];
    end
     if kk == 1
        %plotta il primo segnale forza (non filtrato), la finestra
        %selezionata per filtrare il segnale forza suddetto, il segnale
        %forza filtrato;
%         figure
%         plot (F(:,kk)), hold on, plot(F_filt(:,kk), 'r-.', 'LineWidth', 2), plot(win.*40), xlim([0 400]), grid on,
  
     end
kk=kk+1;
end

end