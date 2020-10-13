
dati= load ('calibrazione.mat')
[g,div_F,div_A,fs] = parametri_fisici();

SNGcal = g * div_A * dati.data(:,2); % Trasformazione del segnale grezzo
time = (1:length (SNGcal))/fs; % creazione del vettore dei tempi

figure (1)
plot (time,SNGcal)
xlabel ('Tempo [s]','FontSize', 18);
ylabel ('Accelerazione [m/s^2]','FontSize', 18);
title  ('Calibrazione', 'FontSize', 18);
%xlim ([0 20])
%ylim ([-20 20])
savefig('Calibrazione.fig') 
% b= timeseries(SNGcal, time);

RMS_SNGcal = rms(SNGcal(round(0.9*end):round(end))); % calcolo del valore RMS

Ref_value = 10; % valore di riferimento per l'RMS dell'accelerazione

CAL = Ref_value / RMS_SNGcal
save('Coefficiente_calibrazione','CAL')

rms (CAL * SNGcal(round(end/3):round(2*end/3))) % RMS del segnale calibrato