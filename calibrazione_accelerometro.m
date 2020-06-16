dati= load ('dati.mat')

[g,div_F,div_A,fs] = parametri_fisici();

SNGcal = g * div_A * dati.data(:,2);

time = (1:length (SNGcal))/fs;

figure (1)
plot (time,SNGcal)
xlabel ('Tempo [s]','FontSize', 18);
ylabel ('Accelerazione [m/s^2]','FontSize', 18);
title  ('Calibrazione', 'FontSize', 18);
xlim ([0 20])
ylim ([-20 20])

b= timeseries(SNGcal, time);