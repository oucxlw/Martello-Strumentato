t=1:2^15;
fs=52100;
t=t/fs;
sigma=0.002;
for i=1:2^15

   sng(i) =  1/(sqrt(2*pi))*exp( -(t(i)-t(2^14))^2 /(2*sigma^2));
end
sng=(sng)';
figure
plot(t,sng)

F=[];
for i=1:30
    F=[F,sng+0.00002*randn(size(sng))];
end


%%
[r,~]=size(sng);
M_win=r;

win=[zeros((r-M_win)/2,1);hann(M_win);zeros((r-M_win)/2,1)];
E_win=sum(win.^2)/length(win);

fft_sng=fft(sng.*win);

fft_sng = (fft_sng)./length (FFT_F(:,1))/fs/E_win; %PSD di F stimata tramite fft
PSD_F_fft(2:end,:)= 2*PSD_F_fft(2:end,:);   