function [ freq, spettro, fase] = smartfftrep( sgn, rep, fs )
%Agisce come smartfft ma mediando rispetto al numero di ripetizioni rep
%indicate

L=length(sgn)/rep;
[freq,  spettro, fase] = smartfft(sgn(1:L), fs);
for i=1:rep-1
    
    [~,  spettrop, fasep] = smartfft(sgn(1+L*i:(i+1)*L), fs);
        spettro=spettro+spettrop;
        fase=fase+fasep;
    end

spettro=spettro/rep;
fase=fase/rep;


end

