function [picchi,n_picchi] = trovacolpi (Forze, soglia, delay, inizio, fine)
% Function che trova i picchi del vettore Forze, usando come Input: 
%I1) il segnale Forze; 
%I2) soglia dei picchi (soglia);
%I3) sample da saltare una volta superata la soglia (delay);
%I4) punto di inizio della ricerca dei picchi (inizio);
%I5) punto di fine della ricerca dei picchi (fine);
%restituendo come Output: 
%O1) il numero di picchi selezionati (picchi_sel); 

x2=Forze;
picchi=[];
ii=inizio;
while   ii < fine(1)
    if  abs(x2(ii)) > soglia
        picchi = [picchi; ii];
        ii=ii+delay;
    end
    ii=ii+1;
end
%Calcolo del numero di picchi
n_picchi = length(picchi);
end