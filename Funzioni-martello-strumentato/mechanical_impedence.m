function [tempstruct] = mechanical_impedence(sng, N,fs)
A={};
tempstruct = struct('MI',A);

% Calcolo ai sensi di quando proposto da:
% Pavement stiffness measurements in relation to mechanical impedance
%
% A partire dall'accelerazine si calcola la Velocià nel tempo e si prende
% il massimo. Com il rapporto tra Fmax e Vmax si calcola MI
% sng è un Nx1 struct che comprende tra gli altri campi Forza (F) e
% Accelerazione (A). Ciascuno di questi è una matrice le cui colonne sono
% le registrazioni dei singoli colpi

[r,~] = size(sng(1).F(:,1));
dt = 1/fs; time1 = 1000*(0:dt:r/fs-dt); %calcolo vettore tempo

%integrazione accelerazione in velocità
for i=1:N
    v_t = cumtrapz (sng(i).A)/fs; %Integriamo a(t)
    % A causa del rumore di fondo si verifica un drift lineare che vogliamo
    % rimuovere. Quindi ci andiamo a prendere il valore di v_t dopo la fine del
    % picco, circa a metà take
    x = round(length(v_t)/2); %prendo la x
    y = v_t(x,:); %ricavo la y
    % costruisco un vettore costante proporzionale all inclinazione
    v_drift = (ones(size(v_t)) .* y/x  );
    % Integro e ricostruisco il drift
    v_drift = cumtrapz (v_drift);
    v_t_corr = v_t - v_drift; %correggo v_t sottraendo il drift
    
    colpo = 2; %A titolo di esempio stampiamo la 2° misura
    
    figure (700), hold on, %stmpo A vs V
    yyaxis left
    plot (time1, sng(i).A_filt(:,colpo)),
    yyaxis right
    %plot(v_t(:,1))
    plot (time1, v_t_corr(:,colpo))
    xlim ([15 40])
    
    figure (701), hold on, %stampo F vs F
    yyaxis left
    plot (time1, sng(i).F_filt(:,colpo)),
    yyaxis right
    plot (time1, v_t_corr(:,colpo))
    xlim ([15 40])
    
    tempstruct(i).MI = max(sng(i).F) ./ max(v_t_corr); %calcolo di MI
    
    figure (702), hold on %plot di MI
    plot (10*log10(tempstruct(i).MI),'*')
    grid on
end
end
