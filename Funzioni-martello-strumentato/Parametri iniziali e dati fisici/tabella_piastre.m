function [piastre] = tabella_piastre()
%UNTITLED genera una tabella con tutte le specifiche delle piastre
%utilizzate
%   Non richiede parametri di imput.
piastre_mass =  [0; 0.006;  0.1967; 0.6274; 0.6249; 1.4293; 2.8871; 15];
piastre_h =     [0; 0.0018; 0.008;  0.008;  0.008;  0.024;  0.0475; 0.16];
piastre_d =     [0; 0.026;  2*sqrt(0.056*0.057/pi); 2*sqrt(0.1*0.1/pi); ...
    2*sqrt(0.1*0.1/pi); 0.1; 0.1;2*sqrt(0.25*0.25/pi)];
piastre = table(piastre_mass,piastre_h,piastre_d);
piastre.Properties.RowNames={'mini','piastrina','quadrata_piccola',...
    'quadrata1','quadrata2','pesante1','pesante2','blocco'};
piastre.Properties.VariableNames={'massa','h','d'}
end

