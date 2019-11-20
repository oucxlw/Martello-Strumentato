function [Matrice] = creamatriceaccelerazione (Segnale, pos, L_pre, L_coda, fs)
% Function che trova i picchi del vettore Forze, usando come Input: 
%I1) il segnale Forze; 
%I2) il segnale Accelerazioni;
%I3) soglia dei picchi (soglia);
%I4) sample da saltare una volta superata la soglia (delay);
%I5) punto di inizio della ricerca dei picchi (inizio);
%I6) punto di fine della ricerca dei picchi (fine);
%I7) lunghezza della parte prima del picco (L_pre);
%I8) lunghezza della parte dopo del picco (L_coda);
%I9) se filt_doppi=1 i colpi vengono filtrati eliminando i doppi colpi (filt_doppi);
%I10) la frequenza di campionamento in Hz;
%restituendo come Output: 
%O1) il numero di picchi selezionati (picchi_sel); 
%O2) la matrice delle singole forze (F) sistemate in colonne;
%O3) la matrice delle singole accelerazioni (A) sistemate in colonne. 

y=Segnale;


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Matrice=[];
%y_sel=[]; %y_sel contencono tutti i picchi selezionati in un unico vettore
n_picchi=length(pos);
for j = 1: n_picchi
    in  = round(pos(j) - L_pre);
    out = round(pos(j) + L_coda);

    y_temp = y(in:out);    

    %y_sel = [y_sel; y_temp];
    Matrice = [Matrice, y_temp];   

j = j+1;
end



end