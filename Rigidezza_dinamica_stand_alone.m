%%
set (0,'DefaultFigureWindowStyle','docked')

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% RIGIDEZZA DINAMICA VASQUEZ stand alone
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

clear variables
Inputs=load('Outputs.mat');
f = Inputs.frequenze.f;
f_fft = Inputs.frequenze.f_fft;
sng = Inputs.sng;
win_1 = Inputs.win_1;
PSD = Inputs.PSD;
FFT = Inputs.FFT;
win_A = Inputs.win_A;
result = Inputs.result;
conf = Inputs.conf;
frequenze = Inputs.frequenze;

[piastre] = tabella_piastre ();
[campioni] = tabella_campioni (conf(1).conf,piastre);

N = length(result.indici_colpo.max_A);
fs = 52100;


for i = 1:N
    m(i) = piastre.massa(conf(i).conf.piastra); %massa della piastra in uso
    h(i) = campioni.h(conf(i).conf.campione);
    s(i) = pi*(campioni.d(conf(i).conf.campione)/2)^2;
end

clear Inputs;

% Plotting
ascissamin = 100;          % Frequenza minima da plottare nei grafici
ascissamax = 2000;       % Frequenza massima da plottare nei grafici

colore = [
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

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Figura confronto Accelerazione-Coerenza secondo Vasquez
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure (405), hold on
title ('Frequenza di risonanza (V.F. Vazquez, S.E. Paje, 2012) ','FontSize',16)
f_analisi=[100 1000];
yyaxis left
hold on
massimi = cell(1,N);
%fr_misura = zeros(1,N);
in = find(f>f_analisi(1),1);
out = find(f<f_analisi(2));
out = out(end);

for i = 1:N
    
    temp = PSD(i).A;
    temp2 = PSD(i).F;
    F_RMS = (sum(temp2,1)/2).^(1/2);
    %F_RMS = (abs((sum(PSD(i).F)), 1)).^(1/2);
    temp_massimi=zeros(1,size(temp,2));
    for j=1:size(temp,2)
        %se la forza è compresa nell'intervalo di forza richiesto scrivo il
        %valore della frequenza di picco, altrimenti lo lascio 0
        %if F_RMS(j)<(3.2) && F_RMS(j)>(1.7) 
        temp_massimi(j) = f(temp(:,j)==max(temp(in:out,j)));%cerco la frequanza di risonanza
        
        plot(f,100*(temp(:,j)./max(temp(in:out,j))),'DisplayName',['Sequenza n°',...
        num2str(i)],'LineWidth',2)
        %end
    end
    %massimi(i) = f(temp==max(temp(in:out)));%cerco la frequanza di risonanza
    result.indici_colpo.max_V(i).fd = temp_massimi;
    %pnt_misura.indici.max_A(i).fd = massimi(i); %scrivo il risultato
    
end
ylim ([0 150])
ylabel('Accelerazione [%]','FontSize',18),

yyaxis  right
hold on
for i = 1:N    
    plot(f_fft,abs(sng(i).Cxy),'DisplayName',['Sequenza n°',num2str(i)],...
        'LineWidth',2)
end
ylabel('Coerenza','FontSize',18),
xlim (f_analisi),ylim([0 1]);
S_primo_0 = zeros(1,N);


figure (406) %S' in funzone della forza RMS
for i=1:N
    result.indici_colpo.max_V(i).K0 = (2*pi* result.indici_colpo.max_V(i).fd ).^2 * m(i);
    result.indici_colpo.max_V(i).S_star = (2*pi* result.indici_colpo.max_V(i).fd ).^2 * m(i)/s(i);
    result.indici_colpo.max_V(i).E  = (2*pi* result.indici_colpo.max_V(i).fd ).^2 * m(i) * h(i)/s(i);   
    
    temp2 = PSD(i).F;
    F_RMS = (sum(temp2,1)/2).^(1/2);
    %plot(max(sng(i).F),  result.indici_colpo.max_V(i).S_star/10^6,'+','LineWidth',2,...
    plot(F_RMS(result.indici_colpo.max_V(i).fd>0),...
        result.indici_colpo.max_V(i).S_star(result.indici_colpo.max_V(i).fd>0)/10^6,...
        '+','LineWidth',2,'color', string(colore(i,:)),'DisplayName',...
        ['Sequenza n°', num2str(i)])

    X = F_RMS(result.indici_colpo.max_V(i).fd>0);
    Y = result.indici_colpo.max_V(i).S_star(result.indici_colpo.max_V(i).fd>0);
    P = polyfit(X,Y,1);
    
    X_exp = 0:max(X*1.5);
    Y_fit = P(1)*X_exp+P(2);
    S_primo_0(i) = Y_fit(1);
    hold on;
    plot(X_exp,Y_fit/10^6,'r-.', 'color', string(colore(i,:)),...
        'DisplayName',['Fit lineare n°', num2str(i)]);
end
grid on;
ylabel('Rigidezza Dinamica apparente x10^6 [N/m]','FontSize',18),
xlabel ('Forza [N]', 'FontSize', 18)
title ('Rigidezza dinamica in funzione della forza impressa','FontSize',20)
xlim=([0 max(F_RMS)])
save('S_primo_0','S_primo_0');
exportgraphics(gcf,'Dstiff_vs_F-Vasquez.pdf', 'BackgroundColor', 'none') 
saveas(gcf, cell2mat(['Dstiff_vs_F-Vasquez_',conf(i).conf.campione,'_',conf(i).conf.piastra,'.fig']))
%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure (402), hold on
title ('Frequenza di risonanza (V.F. Vazquez, S.E. Paje, 2012) ','FontSize',16)

yyaxis left
hold on
fr_misura = zeros(1,N);
for i = 1:N
    temp = abs(PSD(i).Aav);
    in = find(f>f_analisi(1),1);
    out = find(f<f_analisi(2));
    out = out(end);
    fr_misura(i) = f(temp==max(temp(in:out)));%cerco la frequanza di risonanza
    pnt_misura.indici.max_A(i).fd=fr_misura(i); %scrivo il risultato
    plot(f,100*(temp/max(temp(in:out))),'DisplayName',['Sequenza n°',...
        num2str(i),' f_r=',num2str(fr_misura(i)),' Hz'],'LineWidth',2)
end
ylim ([0 150])
ylabel('Accelerazione [%]','FontSize',18),

yyaxis  right
hold on
for i = 1:N    
    plot(f_fft,abs(sng(i).Cxy),'DisplayName',['Sequenza n°',num2str(i)],...
        'LineWidth',2)
end
ylabel('Coerenza','FontSize',18),

xlim (f_analisi),ylim([0 1]);
xlabel('Frequenza [Hz]','FontSize',18),
set(gca, 'XScale', 'log'), %set(gca, 'YScale', 'log'),
legend('FontSize',12)
% Requires R2020a or later
exportgraphics(gcf,'Acc-Coher-Vasquez.pdf', 'BackgroundColor', 'none') 
saveas(gcf, cell2mat(['Acc-Coher-Vasquez_',conf(i).conf.campione,'_',conf(i).conf.piastra,'.fig']))

for i=1:N
    %pnt_misura.indici.max_A(i).K0 = fr_misura(i);
    pnt_misura.indici.max_A(i).K0 = (2*pi*fr_misura(i))^2 * m(i);
    pnt_misura.indici.max_A(i).S_star = (2*pi*fr_misura(i))^2 * m(i)/s(i);
    pnt_misura.indici.max_A(i).E  = (2*pi*fr_misura(i))^2 * m(i) * h(i)/s(i);   
end
s_pnt_misura = vertcat(pnt_misura.indici.max_A.S_star);

save('S_primo','s_pnt_misura');
save('Outputs.mat', 'sng', 'PSD', 'FFT', 'result','frequenze', 'win_1', 'win_A', 'conf');


