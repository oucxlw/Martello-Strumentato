function [radicale,phi,angolo,w] = radice(m,k,c,freq)
%Parte funzionale di A(w) e phi, calcolato anche come angolo

w=2*pi*freq;
w0=sqrt(k/m);

radicale=( (1-(w/w0).^2).^2 + (w.*c/k).^2 ).^0.5;
argatan=( w.*c/k )./( 1-(w/w0).^2 );
phi=-atan(argatan);

angolo=-angle(-w.^2*m+1i*w.*c+k);


end

