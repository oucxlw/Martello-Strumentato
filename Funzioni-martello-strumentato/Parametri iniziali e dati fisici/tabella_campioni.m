function [campioni] = tabella_campioni(conf,piastre)
%UNTITLED2 Tabella specifiche campioni
%   Dare come parametri di ingresso la tabella configurazione conf e la
%   tablla piastre.

campioni_mass = [0.658; 0.634; 0.639; 0.607; 0.541; 0.586; 1; 1; 0.5531; 0.3926; 0.1461; 0.123; 0.6128; 1.446; 0.0383;  0.0705; 0.1064; 1; 1; 1; 1; 1; 1; 1];
campioni_h =    [0.037; 0.035; 0.036; 0.032; 0.029; 0.032; 0.020; 1; 0.031; 0.027;  0.031;  0.031; 0.039;  0.09;  0.005;    0.01;  0.015; 0.045; 0.045; 0.045; 0.019; 0.05; 0.04; 0.04];
campioni_d =    [0.098; 0.098; 0.098; 0.098; 0.098; 0.098; 2*sqrt(0.50*0.24/pi); 1; 0.1;   0.099;  0.097;  0.097;   0.1;   0.1;  2*sqrt(0.098*0.096/pi);...
    0.1; 0.1; 2*sqrt(0.3*0.9/pi); piastre.d(conf.piastra); piastre.d(conf.piastra); ...
    piastre.d(conf.piastra); piastre.d(conf.piastra); piastre.d(conf.piastra); piastre.d(conf.piastra)];
campioni = table(campioni_mass,campioni_h,campioni_d);
campioni.Properties.RowNames={'8_2%','10_2%','11_2%','13_0%','15_0%','16_0%','slab2','calibrazione','c0','c1', 'c2', 'c2b', 'c3','c4','polipropilene','teflon','PVC',...
    'slab','viabattelli','viacocchi','legno','arezzo1','massarosa1','massarosa_b'};
campioni.Properties.VariableNames = {'massa','h','d'} %#ok<NOPRT>
end

