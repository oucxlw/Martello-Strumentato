
load('misura1.mat')
Cxy1=Cxy;
PSD_K_1=PSD_Kav_misura;
devst1=devst_K_sintetico;
FFT_K_1=FFT_K_misura;

load('misura2.mat')
Cxy2=Cxy;
PSD_K_2=PSD_Kav_misura;
devst2=devst_K_sintetico;
FFT_K_2=FFT_K_misura;

load('misura3.mat')
Cxy3=Cxy;
PSD_K_3=PSD_Kav_misura;
devst3=devst_K_sintetico;
FFT_K_3=FFT_K_misura;

%quarta misura opzionale
load('misura4.mat')
Cxy4=Cxy;
PSD_K_4=PSD_Kav_misura;
devst4=devst_K_sintetico;
FFT_K_4=FFT_K_misura;

load('misura5.mat')
Cxy5=Cxy;
PSD_K_5=PSD_Kav_misura;
devst5=devst_K_sintetico;
FFT_K_5=FFT_K_misura;

colore=[
    '.000, .447, .741';
    '.850, .325, .098';
    '.929, .694, .125';
    '.494, .184, .556';
    '.466, .674, .188';
    '.301, .745, .933';
    '.635, .078, .184';
    '.762, .762, .762';
    '.999, .682, .788';
    '.000, .502, .502';
    '.000, .000, .627';
    '.430, .135, .078';
    '.559, .354, .075';
    '.424, .124, .526';
    '.446, .644, .148';
    '.301, .705, .903';
    '.615, .018, .114'];

figure
hold on
subplot(4,1,1), hold on
title('Coerenza')
plot(f_fft,mean(Cxy1,2),'color',string(colore(1,:)),'LineWidth', 2)
plot(f_fft,mean(Cxy2,2),'color',string(colore(2,:)),'LineWidth', 2)
plot(f_fft,mean(Cxy3,2),'color',string(colore(3,:)),'LineWidth', 2)
plot(f_fft,mean(Cxy4,2),'color',string(colore(4,:)),'LineWidth', 2)
plot(f_fft,mean(Cxy5,2),'color',string(colore(5,:)),'LineWidth', 2)

legend('10 N','20 N','100 N','200 N','1000')

xlim([20 5000]),grid on
set(gca, 'XScale', 'log'),

subplot(4,1,[2,3]), hold on
title('Rigidezza dinamica - Modulo')

plot(f,10*log10(mean(abs(PSD_K_1),2)),'color',string(colore(1,:)),'LineWidth', 2)
plot(f,10*log10(mean(abs(PSD_K_2),2)),'color',string(colore(2,:)),'LineWidth', 2)
plot(f,10*log10(mean(abs(PSD_K_3),2)),'color',string(colore(3,:)),'LineWidth', 2)
plot(f,10*log10(mean(abs(PSD_K_4),2)),'color',string(colore(4,:)),'LineWidth', 2)
plot(f,10*log10(mean(abs(PSD_K_5),2)),'color',string(colore(5,:)),'LineWidth', 2)

plot(f,10*log10(mean(abs(PSD_K_1),2)+devst1),'--','color',string(colore(1,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_2),2)+devst2),'--','color',string(colore(2,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_3),2)+devst3),'--','color',string(colore(3,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_4),2)+devst4),'--','color',string(colore(4,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_5),2)+devst5),'--','color',string(colore(5,:)),'LineWidth', 1)

plot(f,10*log10(mean(abs(PSD_K_1),2)-devst1),'--','color',string(colore(1,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_2),2)-devst2),'--','color',string(colore(2,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_3),2)-devst3),'--','color',string(colore(3,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_4),2)-devst4),'--','color',string(colore(4,:)),'LineWidth', 1)
plot(f,10*log10(mean(abs(PSD_K_5),2)-devst5),'--','color',string(colore(5,:)),'LineWidth', 1)

xlim([20 5000]),grid on
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),

subplot(4,1,4), hold on
title('Rigidezza dinamica - Fase')

plot(f_fft, 180/pi*unwrap(angle(mean(FFT_K_1,2))),'color',string(colore(1,:)),'LineWidth', 2)
plot(f_fft, 180/pi*unwrap(angle(mean(FFT_K_2,2))),'color',string(colore(2,:)),'LineWidth', 2)
plot(f_fft, 180/pi*unwrap(angle(mean(FFT_K_3,2))),'color',string(colore(3,:)),'LineWidth', 2)
plot(f_fft, 180/pi*unwrap(angle(mean(FFT_K_4,2))),'color',string(colore(4,:)),'LineWidth', 2)
plot(f_fft, 180/pi*unwrap(angle(mean(FFT_K_5,2))),'color',string(colore(5,:)),'LineWidth', 2)

xlim([20 5000]),ylim([-200 200]), yticks([-180 -90 0 90 180]),grid on
set(gca, 'XScale', 'log')

CXY=[Cxy1 Cxy2 Cxy3 Cxy4 Cxy5];
for j= 1:5 
    for i=1:3
        f0_temp(i) = f_fft(find(abs(CXY(:,i+j)) > 0.6, 1));
        fend_temp(i)=f_fft(find(abs(CXY(find(f_fft==f0_temp(1)):end,i+j))<0.6,1))
    end
    f0(j) = mean(f0_temp);
    Df(j) = mean(fend_temp-f0_temp);
end
f0
Df