function [soglia,delay,inizio,fine] = parametri_ricerca_picchi(fs)
%UNTITLED4 Parametri per la ricerca dei picchi
soglia=10;              % Soglia dei picchi;
delay=round(0.5*fs);    % Sample da saltare una volta superata la soglia
inizio=1*fs;            % Punto di inizio della ricerca dei picchi;
fine=round(0.95*size(x));           % Punto di fine della ricerca dei picchi

end

