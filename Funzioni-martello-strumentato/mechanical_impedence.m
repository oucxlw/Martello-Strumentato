function [tempstruct,mdls,MI_av,MI_av_dev,MI_lin,MI_lin_R,MI_exp,MI_exp_R] = mechanical_impedence(sng, N,fs)
A={};
tempstruct = struct('MI',A);

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

% Calcolo ai sensi di quando proposto da:
% Pavement stiffness measurements in relation to mechanical impedance
%
% A partire dall'accelerazine si calcola la Velocià nel tempo e si prende
% il massimo. Com il rapporto tra Fmax e Vmax si calcola MI
% sng è un Nx1 struct che comprende tra gli altri campi Forza (F) e
% Accelerazione (A). Ciascuno di questi è una matrice le cui colonne sono
% le registrazioni dei singoli colpi

[r,~] = size (sng(1).F(:,1));
dt = 1/fs; time1 = 1000*(0:dt:r/fs-dt); %calcolo vettore tempo
y_max = 0; %inizializzazione della y max per i limiti del plot
x_max = 0; %inizializzazione della x max per i limiti del plot
X=[]; %vettore delle X per il fit
Y=[]; %vettore delle Y per il fit
MI_av = zeros(1,N);
MI_lin = zeros(1,N);
MI_exp = zeros(1,N);
MI_av_dev = zeros(1,N);
MI_lin_R = zeros(1,N);
MI_exp_R = zeros(1,N);
filtro_intensita = 0; %Serve a togliere dal fit misure con Fmax troppo elevata
lim_Fmax = 700; %Forza a cui si taglia per il fit
%integrazione accelerazione in velocità e plot
for i=1:N
    v_t = cumtrapz (sng(i).A)/fs; %Integriamo a(t)
    % A causa del rumore di fondo si verifica un drift lineare che vogliamo
    % rimuovere. Quindi ci andiamo a prendere il valore di v_t dopo la fine del
    % picco, circa a metà take
    x = round (length (v_t)/2); %prendo la x
    y = v_t (x,:); %ricavo la y
    % costruisco un vettore costante proporzionale all inclinazione
    v_drift = (ones (size (v_t)) .* y/x  );
    % Integro e ricostruisco il drift
    v_drift = cumtrapz (v_drift);
    v_t_corr = v_t - v_drift; %correggo v_t sottraendo il drift
    
    tempstruct(i).MI = max(sng(i).F) ./ max(abs(v_t_corr)); %calcolo di MI
    tempstruct(i).MI = max(sng(i).F) ./ max(v_t_corr); %calcolo di MI

    %Cerco i limiti del grafico 702
    if y_max < max (tempstruct(i).MI)
        y_max = max (tempstruct(i).MI);
    end
    if x_max < length (tempstruct(i).MI)
        x_max = length (tempstruct(i).MI);
    end
    
    colpo = 2; %A titolo di esempio stampiamo la 2° misura
    
    figure (700), hold on, %stmpo A vs V
    yyaxis left
    plot (time1, sng(i).A_filt(:,colpo)),
    yyaxis right
    %plot(v_t(:,1))
    plot (time1, v_t_corr(:,colpo))
    
    figure (701), hold on, %stampo F vs F
    yyaxis left
    plot (time1, sng(i).F_filt(:,colpo)),
    yyaxis right
    plot (time1, v_t_corr(:,colpo))
    
    % Plot degli indicatori MI dei singoli impatti
    figure (702), hold on %plot di MI
    plot (20*log10(tempstruct(i).MI),'*', 'color', string(colore(i,:)),...
        'DisplayName',['Sequenza n°',num2str(i)])
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Plot degli MI in funzione del picco di forza
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    figure (703), hold on %plot di MI
    plot (max(sng(i).F),20*log10(tempstruct(i).MI),'+','LineWidth',2,...
        'color', string(colore(i,:)),'DisplayName', ['Sequenza n°', num2str(i)])
    
    X = max(sng(i).F);
    Y = 20*log10(tempstruct(i).MI);
    
    % filtro intensità
    if filtro_intensita==1
        Y(X>lim_Fmax) = [];
        X(X>lim_Fmax) = [];
    end
    
    % fit esponenziale
    modelfun = @(b,x)b(1)*exp(-(x-b(4))/b(2))+b(3);
    mdls(i).mdl_exp = fitnlm(X, Y, modelfun, [-5 100 mean(Y(round(2/3*end):end)) 0]); % Modello dei primi picchi
    %[fit,R] = nlinfit(X, Y, modelfun, [200 150 150 0]);
    X_exp = 1:max(X*1.1);
    Y_exp = feval(mdls(i).mdl_exp,X_exp);
    % plot della curva esponenziale
    plot (X_exp, feval(mdls(i).mdl_exp,X_exp), 'color', string(colore(i,:)),...
        'LineWidth',2,'DisplayName',['Fit exp n°', num2str(i),' R^2=',num2str(...
        mdls(i).mdl_exp.Rsquared.Adjusted)])% plot del modello
    MI_exp(i) = Y_exp(1);
    MI_exp_R(i) = mdls(i).mdl_exp.Rsquared.Adjusted;
    
    % fit lineare
    modelfun = @(b,x)b(1)*x+b(2);
    mdls(i).mdl_lin = fitnlm(X, Y, modelfun, [1 0]); % Modello dei primi picchi
    P = polyfit(X,Y,1);
    Y_fit = feval(mdls(i).mdl_lin,X_exp);
    hold on;
    % plot della retta
%     plot(X_exp,Y_fit,'r-.', 'color', string(colore(i,:)),...
%         'LineWidth',2,'DisplayName',['Fit lineare n°', num2str(i),' R^2=',num2str(...
%        mdls(i).mdl_lin.Rsquared.Adjusted)]);
    % eguaglio l'impedenza meccanica all'intercetta
    MI_lin(i) = mdls(i).mdl_lin.Coefficients.Estimate(2);
    MI_lin_R(i) = mdls(i).mdl_lin.Rsquared.Adjusted;

    % media
    MI_av(i) = mean(20*log10(tempstruct(i).MI));
    MI_av_dev(i) = std(20*log10(tempstruct(i).MI));
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % plot lineare non dB dell'impedenza meccanica
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    figure (706), hold on
    plot (max(sng(i).F),(tempstruct(i).MI),'+','LineWidth',2,...
        'color', string(colore(i,:)),'DisplayName', ['Sequenza n°', num2str(i)])
    
    X = max(sng(i).F);
    Y2 = tempstruct(i).MI;
    
    % filtro intensità
    if filtro_intensita==1
        Y2(X>lim_Fmax) = [];
        X(X>lim_Fmax) = [];
    end
    
    % fit esponenziale
%     modelfun = @(b,x)b(1)*exp(-(x-b(4))/b(2))+b(3);
%     mdls(i).mdl_exp_bis = fitnlm(X, Y2, modelfun, [-5 100 mean(Y(round (2/3*end):end)) 0]); % Modello dei primi picchi
%     %[fit,R] = nlinfit(X, Y, modelfun, [200 150 150 0]);
%     X_exp = 1:max(X*1.1);
%     Y_exp = feval(mdls(i).mdl_exp_bis,X_exp);
%     % plot della curva esponenziale
%     plot (X_exp, feval(mdls(i).mdl_exp_bis,X_exp), 'color', string(colore(i,:)),...
%         'LineWidth',2,'DisplayName',['Fit exp n°', num2str(i),' R^2=',num2str(...
%         mdls(i).mdl_exp_bis.Rsquared.Adjusted)])% plot del modello
    % MI_exp(i) = Y_exp(1);
    % MI_exp_R(i) = mdls(i).mdl_exp_bis.Rsquared.Adjusted;
    % fit lineare
    modelfun = @(b,x)b(1)*x+b(2);
    mdls(i).mdl_lin_bis = fitnlm(X, Y2, modelfun, [1 0]); % Modello dei primi picchi
    P = polyfit(X,Y,1);
    Y_fit = feval(mdls(i).mdl_lin_bis,X_exp);
    hold on;
    % plot della retta
    plot(X_exp,Y_fit,'r-.', 'color', string(colore(i,:)),...
        'LineWidth',2,'DisplayName',['Fit lineare n°', num2str(i),' R^2=',num2str(...
        mdls(i).mdl_lin_bis.Rsquared.Adjusted)]);
    % eguaglio l'impedenza meccanica all'intercetta
    %MI_lin(i) = mdls(i).mdl_lin.Coefficients.Estimate(2);
    %MI_lin_R(i) = mdls(i).mdl_lin.Rsquared.Adjusted;
end

figure (700)
yyaxis left
ylabel ('Accelerazione','FontSize', 18)
yyaxis right
ylabel ('Velocità','FontSize', 18)
xlim ([15 40])

figure (701)
yyaxis left
ylabel ('Forza','FontSize', 18)
yyaxis right
ylabel ('Velocità','FontSize', 18)
xlim ([15 40])

figure (702), hold on %plot di MI
grid on
grid minor
xlim ([0 x_max+1])
xticks ( 1:x_max);
ylim ([0 1.2*10*log10(y_max)])
xlabel ('Colpo progressivo', 'FontSize', 18)
ylabel ('Impedenza meccanica MI [dB @ 1Ns/m]', 'FontSize', 18)
legend ('FontSize', 18,'Location','southeast')

figure (703), hold on %plot di MI vs max F
title('Impedenza Meccanica vs Massimo della forza','FontSize',20)
grid on
grid minor
%xlim ([0 x_max+1])
%xticks ( 1:x_max);
%ylim ([0 1.2*10*log10(y_max)])
xlabel ('Forza [N]', 'FontSize', 18)
ylabel ('Impedenza meccanica MI [dB @ 1Ns/m]', 'FontSize', 18)
legend ('FontSize', 12,'Location','northwest')
saveas (gcf, 'MI_vs_F_M_Li_2016.fig')
exportgraphics (gcf,'MI_vs_F_M_Li_2016.pdf', 'BackgroundColor', 'none') 

figure (704), hold on
bar (1, MI_exp);%,'DisplayName', ['Fit esponenziale n°', num2str(i)]);
%bar (1, MI_av_lin);%,'DisplayName', ['Fit esponenziale n°', num2str(i)]);
title('Impedenza Meccanica ottenuta con fit esponenziale','FontSize',20)
xticks ([])
ylabel ('Impedenza meccanica MI'' [dB @ 1Ns/m]', 'FontSize', 18)
grid on
grid minor
legend ('FontSize', 12)
saveas (gcf, 'Impedenza_meccanica-M_Li_2016.fig')
exportgraphics (gcf,'Impedenza_meccanica-M_Li_2016.pdf', 'BackgroundColor', 'none') 

figure (705), hold on
%bar (1, MI_av);%,'DisplayName', ['Fit esponenziale n°', num2str(i)]);
bar (1, MI_lin);%,'DisplayName', ['Fit esponenziale n°', num2str(i)]);
title('Impedenza Meccanica ottenuta con fit lineare','FontSize',20)
xticks ([])
ylabel ('Impedenza meccanica MI'' [dB @ 1Ns/m]', 'FontSize', 18)
grid on
grid minor
legend ('FontSize', 12)
saveas (gcf, 'Impedenza_meccanica-M_Li_2016_lin.fig')
exportgraphics (gcf,'Impedenza_meccanica-M_Li_2016_lin.pdf', 'BackgroundColor', 'none')

figure (706), hold on %plot non decibel di MI vs max F
title('Impedenza Meccanica vs Massimo della forza','FontSize',20)
grid on
grid minor
%xlim ([0 x_max+1])
%xticks ( 1:x_max);
%ylim ([0 1.2*10*log10(y_max)])
xlabel ('Forza [N]', 'FontSize', 18)
ylabel ('Impedenza meccanica MI [Ns/m]', 'FontSize', 18)
legend ('FontSize', 12,'Location','southeast')
saveas (gcf, 'MI_vs_F_M_Li_2016_linear.fig')
exportgraphics (gcf,'MI_vs_F_M_Li_2016_linear.pdf', 'BackgroundColor', 'none') 
end