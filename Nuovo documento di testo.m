% Martello strumentato filtro intensita
set (0,'DefaultFigureWindowStyle','docked')
clc
close all
clear variables
load dati.mat

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Importazione di forza e accelerazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
campione='teflon senza biadesivo, rip 1';
accelerometro=0;
punta='m'; %m=metallica; p=plastica; g=gomma.
piastra='pesante';
martellamento='auto';

x = pp_m_teflon_1 (:,1); % Force [N] 
y = pp_m_teflon_1 (:,accelerometro+2); % Accelerazione [m/s^2] 

%x = polipropilene_piccola_m_biad_3 (:,1); % Force [N] 
%y = polipropilene_piccola_m_biad_3 (:,accelerometro+2); % Accelerazione [m/s^2] 

%  x = slab_piccola_m_res_4 (:,1); % Force [N] 
%  y = slab_piccola_m_res_4 (:,accelerometro+2); % Accelerazione [m/s^2] 
%x = sito_piastra_piccola_lato_metallo_mart; % Force [N] 
%y = sito_piastra_piccola_lato_metallo_acc; % Accelerazione [m/s^2] 

