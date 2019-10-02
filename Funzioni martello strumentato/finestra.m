function [F_filt, A_filt, V_filt, D_filt, F_filt_all, A_filt_all] = finestra (Forze, Accelerazioni, picchi_sel, window_F, window_A, fs)
%Function che esegue la finestrature delle colonne delle matrici Forza e Accelerazioni,
%permettendo di scegliere diverse tra varie finestre.
%Prende in Input: 
%I1) la matrice delle Forze (F); ù
%I2) la matrice delle accelerazioni (A); 
%I3) il numero di picchi selezionati (picchi_sel);
%I4) il tipo di finestra selezionata per le Forze (window_F=numero tra 0 e 5), scegliendo tra: 
        %case 0 = nessuna finestratura,
        %case 1 = onda quadra come finestra,
        %case 2 = funzione di hamming come finestra,
        %case 3 = funzione di Blackman-Harris (versione 1) come finestra,
        %case 4 = funzione di Blackman-Harris (versione 2) come finestra,
        %case 5 = funzione di hann come finestra;
%I5) il tipo di finestra selezionata per le Accelerazioni (window_A=numero tra 0 e 3), scegliendo tra: 
        %case 0 = nessuna finestratura,
        %case 1 = funzione di esponenziale come finestra,
        %case 2 = funzione di Blackman-Harris come finestra,
        %case 3 = funzione di hann come finestra;
%Restituendo come Output: 
%O1) la matrice Forze dopo il finestratura (F_filt);
%O2) la matrice Accelerazioni dopo il finestratura (A_filt);
%O3) la matrice Velocità (derivate da A_filt) dopo il finestratura (V_filt);
%O4) la matrice Spostamenti (derivate da V_filt) dopo il finestratura (D_filt);
%O5) il vettore Forze dopo il finestratura (F_filt_all) necessario per il calcolo della funzione coerenza;
%O6) il vettore Accelerazioni dopo il finestratura (A_filt_all) necessario per il calcolo della funzione coerenza.


F=Forze;
A=Accelerazioni;

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione???
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
    switch window_F
        case 0
            F_filt = [F_filt,F(:, kk)];
        
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
            Sgn=Sgn.*win;
            zeri=zeros(length(F(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            F_filt=[F_filt, temp];
    end
    
% % %     %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% % %     % Calcolo NORMA F per il singolo colpo
% % %     %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% % %     if norm==1
% % %         temp=0;
% % %         for ll=1:length(F(:,kk))
% % %             temp=temp+(F(ll,kk)^2);        
% % %         end
% % %         norma_F = temp^(1/2);
% % %     else
% % %         norma_F = 1;
% % %     end
% % %     F_filt(:,kk) = F_filt(:,kk)./norma_F(end);
% % %     F_filt_all = [F_filt_all; F_filt(:,kk)]; 

    %<<<<<<<<<<<<<<<<<<<
    % Finestratura di A
    %<<<<<<<<<<<<<<<<<<<
    switch window_A
        case 0 %nessuna finestra
            A_filt= [A_filt, A(:, kk)];      
            A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A(:,kk)];
            
            V_filt=[V_filt, cumtrapz(time,A(:,kk))];       %Calcolo Velocità V a partire da A non finestrato
            D_filt=[D_filt, cumtrapz(time,V_filt(:, kk))];    %Calcolo Spostamento D a partire da V non finestrato
% % %             V_filt(:,kk) = V_filt(:,kk)/norma_F(end);       %NON finestro e nel caso normalizzo V
% % %             D_filt(:,kk) = D_filt(:,kk)/norma_F(end);       %NON finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];          %Accodo D
            
        case 1% finestra: Esponenziale:           
            Lw2=length(A(:,1));
            tw2=0:1:Lw2-1;
            beta=50; %<<<<<<<<<<<<<<<<<
            w2=1.5*exp(-tw2./beta);
            y_sel_filt = A(:, kk).*w2';
            A_filt= [A_filt, y_sel_filt];                 
            A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A_filt(:,kk)];      
       
        case 2 % finestra: Blackman-Harris
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
            D_filt=[D_filt, cumtrap((1/fs),V_filt(:, kk))];    %Calcolo Spostamento D a partire da V non finestrato
% % %             V_filt(:,kk) = win.*V_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo V
% % %             D_filt(:,kk) = win.*D_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];          %Accodo D
            
        case 3 % finestra: hann  
            t0=1; %L_pre-3; % inizio della finestra di salita
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
% % %             A_filt(:,kk) = A_filt(:,kk)/norma_F(end);
            A_filt_all=[A_filt_all; A_filt(:,kk)];      
            Lung=length(A_filt(:,1));
            f =0:fs/Lung:fs/2;

            V_filt=[V_filt, cumtrapz(time,A(:,kk))];       %Calcolo Velocità V a partire da A non finestrato
            D_filt=[D_filt, cumtrapz(time,V_filt(:,kk))];    %Calcolo Spostamento D a partire da V non finestrato
% % %             V_filt(:,kk) = win.*V_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo V
% % %             D_filt(:,kk) = win.*D_filt(:,kk)/norma_F(end);  %finestro e nel caso normalizzo D
            V_filt_all=[V_filt_all; V_filt(:,kk)];          %Accodo V
            D_filt_all=[D_filt_all; D_filt(:,kk)];          %Accodo D      
    end
    
    if kk == 1
        %plotta il primo segnale forza (non filtrato), la finestra
        %selezionata per filtrare il segnale forza suddetto, il segnale
        %forza filtrato;
        %come sopra, per l'accelerazione.
        figure
        subplot(2,1,1), plot (F(:,kk)), hold on, plot(F_filt(:,kk), 'r-.', 'LineWidth', 2), plot(win.*40), xlim([0 400]), grid on,
        subplot(2,1,2), plot(A(:,kk)), hold on, plot(A_filt(:,kk), 'r-.', 'LineWidth', 2), plot(win.*100),xlim([0 1000]), grid on,
    end
kk=kk+1;
end

end