%%
% Confronti tra ripetizioni
%clear variables
load dati_confronto_materiali.mat
close all

figure (1)
hold on 
titolo= 'Confronto materiali, lato';
sgtitle(titolo)
% colori: #0072BD #D95319 #EDB120 #7E2F8E #77AC30 #7E2F8E

%confronto campioni di materiali diversi
semilogx (f, 20*log10(Dstiff_av_teflon),        'LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_PVC),           'LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_polipropilene), 'LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_C3),            'LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_teflonPVC),     'LineWidth', 2),

semilogx (f, 20*log10(Dstiff_av_teflon_6_10),        '--','color','#0042BD','LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_PVC_6_10),           '--','color','#D95319','LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_polipropilene_6_10), '--','color','#EDB120','LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_C3_6_10),            '--','color','#7E2F8E','LineWidth', 2),
semilogx (f, 20*log10(Dstiff_av_teflonPVC_6_10),     '--','color','#77AC30','LineWidth', 2),

% semilogx (f, 20*log10(Dstiff_av_teflon_1_5),        ':','color','#0072BD','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff_av_PVC_1_5),           ':','color','#D95319','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff_av_polipropilene_1_5), ':','color','#EDB120','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff_av_C3_1_5),            ':','color','#7E2F8E','LineWidth', 2),
% semilogx (f, 20*log10(Dstiff_av_teflonPVC_1_5),     ':','color','#77AC30','LineWidth', 2),


legend('Teflon 1 cm','PVC 1.5 mm','Polipropilene 0.5 mm','Campione 3','Teflon + PVC',...
    'Teflon 1 cm (forti)','PVC 1.5 mm (forti)','Polipropilene 0.5 mm (forti)','Campione 3 (forti)',...
    'Teflon + PVC (forti)')


set(gca, 'Xscale', 'log'), set(gca, 'YScale', 'log'),
%ylim([100 200])
xlabel('log(Frequency) [Hz]')
ylabel('20 log |Dynamic Stiffness| (dB ref 1 N/m]'),
%title(['Dynamic Stiffness (Force/Displacement) Amplitude']), 
grid on,
xlim([20 2000]), %ylim([130 170])
hold off
%legend('Piccola biad','Piccola Attak','Piccola pattex','Grande biadesivo','Grande attak','Grande pattex','Grande pattex 2gg')

%legend('Grande Rip. 1','Grande Rip. 2','Grande Rip. 3','Piccola Rip. 1','Piccola Rip. 2','Piccola Rip. 3','Grande Rip. 1 no av','Grande Rip. 2 no av','Grande Rip. 3 no av','Piccola Rip. 1 no av','Piccola Rip. 2 no av','Piccola Rip. 3 no av')
saveas (gcf, [titolo,'.fig'])




K1= 116; K2= 129; K3= 140; K4= 141; K5= 122;
h1 = 0.01; h2= 0.015; h3= 0.005; h4= 0.01; h5= 0.025;
H=[h1 h2 h3 h4 h5]
r1 = 0.0025; r2= 0.0025; r3= 0.0025; r4= 0.0025; r5= 0.0025;

s1= pi*r1^2; s2= pi*r2^2; s3= pi*r3^2; s4= pi*r4^2; s5= pi*r5^2;
S=[s1 s2 s3 s4 s5]
E_teflonPTFE= (10^(K1/20)*h1)/(s1)/10^9
E_PVC       = (10^(K2/20)*h2)/(s1)/10^9
E_poliprop  = (10^(K3/20)*h3)/(s1)/10^9
E_C3        = (10^(K4/20)*h4)/(s1)/10^9
E5          = (10^(K5/20)*h5)/(s1)/10^9

E_teflonPTFE_e= 0.4
E_PVC_e       = (2.4)
E_poliprop_e  = (1.5)
E_C3_e        = 0.5
E5_e          = 1

E_teflonPTFE_e2= 0.4
E_PVC_e2       = 4.1
E_poliprop_e2  = 2
E_C3_e2        = 1.5
E5_e2          = 1

x=categorical( {'Teflon', 'PVC','Polipropilene', 'C3','Teflon + polipropilene' })

x=reordercats(x,{'Teflon', 'PVC','Polipropilene', 'C3','Teflon + polipropilene' })

Eo = [E_teflonPTFE, E_PVC, E_poliprop, E_C3, E5]
Ee = [E_teflonPTFE_e, E_PVC_e, E_poliprop_e, E_C3_e, E5_e]
Ee2 = [E_teflonPTFE_e2, E_PVC_e2, E_poliprop_e2, E_C3_e2, E5_e2]
EE(1,:)=Ee,EE(2,:)=Ee2,


figure (2), hold on
%bar (x,Eo,0.4)


f1=80;f2=150;
x1=find(f>f1);x1=x1(1);
x2=find(f<f2);x2=x2(end);

DS1 = Dstiff_av_teflon_1_5(x1:x2);
DS2 = Dstiff_av_PVC_1_5(x1:x2);
DS3 = Dstiff_av_polipropilene_1_5(x1:x2);
DS4 = Dstiff_av_C3_1_5(x1:x2);
DS5 = Dstiff_av_teflonPVC_1_5(x1:x2);

data = [DS1, DS2, DS3, DS4, DS5];
for k=1:5
moduli(:,k) = (data(:,k)*H(k))/(s1)/10^9
end


boxplot(EE,'colors',[0 0.5 0])
boxplot(moduli)
ylim=([0 17])


(0.67*(2*pi*57.97)^2)/(pi*0.05^2)
