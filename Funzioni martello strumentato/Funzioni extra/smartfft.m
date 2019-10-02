function [freq,  spettro, fase] = smartfft(sgn, fs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
freq = 0:fs/length(sgn):fs/2;
spettro = 2*abs(fft(sgn))/length(sgn);
spettro = spettro(1:length(freq));
fase=complexphase(fft(sgn));
fase=fase(1:length(freq));
end

