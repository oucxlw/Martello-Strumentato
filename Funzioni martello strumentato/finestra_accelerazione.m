function [A_filt] = finestra_accelerazione (Accelerazione, window_A, L_win, fs)
%Function che esegue la finestratura delle colonne della matrice Accelerazione,
%permettendo di scegliere diverse tra varie finestre.
%Prende in Input: 
%I1) la matrice delle accelerazioni (A); 
%I2) il tipo di finestra selezionata per le Accelerazioni (window_A=numero tra 1 e 3), scegliendo tra: 
        %case 1 = funzione di esponenziale come finestra,
        %case 2 = funzione di Blackman-Harris come finestra,
        %case 3 = funzione di hann come finestra;
        %case 4= finestra come case 5 della forza lunga 2/3 L_win;
        %case 5= finestra come case 5 della forza lunga 2 volte L_win;
%I3) la lunghezza del filtro usato nella forza;  
%I4)la frequenza di campionamento in Hz;
%Restituendo come Output: 
%O1) la matrice Accelerazioni dopo il finestratura (A_filt);


A=Accelerazione;

%<<<<<<<<<<<<<<<
% Finestratura 
%<<<<<<<<<<<<<<<
A_filt=[];% Matrici Forza e Accelerazione picchi selezionati
L = length(A(:,1));
dt=1/fs; time=1000*(0:dt:L/fs-dt);
[r,picchi_sel]=size(A);
for kk=1:picchi_sel
  
    %sottrazione della media
    A(:,kk)=A(:,kk)-mean(A(:,kk));

    %<<<<<<<<<<<<<<<<<<<
    % Finestratura di A
    %<<<<<<<<<<<<<<<<<<<
    switch window_A 
        case 1% finestra: Esponenziale:           
            Lw2=length(A(:,1));
            tw2=0:1:L_win(kk);
            beta=20; %<<<<<<<<<<<<<<<<<
            win=(4*exp(-tw2./beta))';
            zeri=zeros(Lw2-L_win(kk)-1,1);
            win=[win; zeri];
            y_sel_filt = A(:, kk).*win;
            A_filt= [A_filt, y_sel_filt];
       
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
       
        case 4 % finestra: hann basata su hann della forza lunga 2/3 L_win
            t0=1; 
            [~,pos] = max(A(:,kk));
            tf=round(L_win(kk)*2/3);
            Lp = tf - t0;
            Delta_f = round(Lp/5);
            Sgn=A(1:tf-1,kk);
            win0= hann(2*t0);
            win0=win0(1:t0);
            uno = ones(Lp-2-Delta_f,1);
            winf= hann (Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn=Sgn.*win;
            zeri=zeros(length(A(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            A_filt=[A_filt, temp];
             
        case 5 % finestra: hann lunga 2 volte L_win
            t0=1; 
            [~,pos] = max(A(:,kk));
            tf=round(3*L_win(kk));
            Lp = tf - t0;
            Delta_f = round(Lp/5);
            Sgn=A(1:tf-1,kk);
            win0= hann(2*t0);
            win0=win0(1:t0);
            uno = ones(Lp-2-Delta_f,1);
            winf= hann (Delta_f*2);
            winf=winf(Delta_f:end);
            win=[ win0;uno; winf];
            Sgn=Sgn.*win;
            zeri=zeros(length(A(:,kk))-length(Sgn),1);
            temp=[Sgn;zeri];
            A_filt=[A_filt, temp];
    end
    
%     if kk == 1
%         %plotta il primo segnale accelerazione (non filtrato), la finestra
%         %selezionata per filtrare il segnale accelerazione suddetto, il segnale
%         %filtrato;
%         figure
%         plot(A(:,kk)), hold on, plot(A_filt(:,kk), 'r-.', 'LineWidth', 2), plot(win.*100),xlim([0 1000]), grid on,
%     end
kk=kk+1;
end

end