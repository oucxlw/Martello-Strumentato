fs = 52100;
T=1;
t=(1/fs):(1/fs):T; 
w0=1000*2*pi;
sng=zeros(size(t));
%sng=sin(w0*t)+1./2*sin(2*w0*t)+1./3*sin(8*w0*t)+cos(1.5*w0*t);%+(1+cos(w0*t/2)).*randn(size(t));
for i=1:15
    
    sng=sng+sin(i*fs/10*t);%*(1/50^3)*(y-50)^3;
end
close all
Y=2*fft(sng)/length(t);
f=0:fs-1;

figure(1),hold on,
plot(f(1:end/2+1),abs(Y(1:end/2+1)))
% figure, plot(t,sng)

Pxx=abs(Y).^2;
plot(f(1:end/2+1),Pxx(1:end/2+1))

rect_win=ones(size(sng));
[PSDxx,f_psd]=periodogram(sng, rect_win, length(t),fs,'power')

plot (f_psd,PSDxx,'b*')

Cxx=2*xcorr(sng)/length(t);
Pxxcorr=2*fft(Cxx)/length(t);
semilogx((abs(Pxxcorr)));ylim([0 1])


legend('abs(Y) = abs(2*fft(sng)/length(sng))','abs(Y)^2','periodogram(sng, rect_ win, length(t),fs,"power")')