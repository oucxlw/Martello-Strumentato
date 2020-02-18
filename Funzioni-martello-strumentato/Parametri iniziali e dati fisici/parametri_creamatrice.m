function [Lsample,L_pre,L_coda] = parametri_creamatrice()
%parametri_creamatrice Dimensioni dei campioni e disposizione degli stessi
%nel take.
Lsample=2^11;
L_pre=round((Lsample/16)); % Lunghezza della parte prima del picco
L_coda=round(Lsample-L_pre-1);     % Lunghezza della coda dei segnali
end

