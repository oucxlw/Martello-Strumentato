%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Impedenza meccanica stand alone
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%<<<<<<<<<<<<<<
% intestazione
%<<<<<<<<<<<<<<
set (0,'DefaultFigureWindowStyle','docked')
clear variables
close all
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

%<<<<<<<<<<<<<<<<<<<<<<
% Mechanical impedence
%<<<<<<<<<<<<<<<<<<<<<<

% Calcolo dell'impedenza meccanica come suggerito da
% Pavement stiffness measurements in relation to mechanical impedance
% Mingliang Li, Wim van Keulen, Halil Ceylan, Dongwei Cao,
% mechanical_impedence plotta anche grafici di velocità vs accelerazione
% e velocità vs forza

[tempstruct,MI_mdls,MI_av,MI_av_dev,MI_lin,MI_lin_R,MI_exp,MI_exp_R] = mechanical_impedence(sng, N,fs);

MI_mean = zeros(1,N);
% Applico a result le MI calcolate in mechanical_impedence
for i = 1:N
    result.indici_colpo.max_V(i).MI = tempstruct(i).MI;
    MI_mean(i) = mean(result.indici_colpo.max_V(i).MI);
end
clear tempstruct
F_max_mean = zeros(1,N);
for i=1:N
    result.indici_colpo.max_V(i).MI_av = MI_av(i); %mean(result.indici_colpo.max_V(i).MI);
    result.indici_colpo.max_V(i).MI_lin = MI_lin(i);
    result.indici_colpo.max_V(i).MI_exp = MI_exp(i);
    F_max_mean(i) = mean(max(sng(i).F_filt));
end

save('F_max_mean.mat','F_max_mean');

save('MI_modelli_fit.mat','MI_mdls');
save('MI_mean.mat','MI_av');
save('MI_exp.mat','MI_exp');
save('MI_lin.mat','MI_lin');

save('Outputs.mat', 'sng', 'PSD', 'FFT', 'result','frequenze', 'win_1', 'win_A', 'conf');


