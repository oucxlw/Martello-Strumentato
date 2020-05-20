set (0,'DefaultFigureWindowStyle','docked')

clear variables
Inputs  =load('Outputs.mat');
f       = Inputs.frequenze.f;
f_fft   = Inputs.frequenze.f_fft;
sng     = Inputs.sng;
win_1   = Inputs.win_1;
PSD     = Inputs.PSD;
FFT     = Inputs.FFT;
win_A   = Inputs.win_A;
result  = Inputs.result;
conf    = Inputs.conf;

[piastre]   = tabella_piastre ();
[campioni]  = tabella_campioni (conf,piastre);

m = piastre.massa(conf.piastra); %massa della piastra in uso
h = campioni.h(conf.campione);
s = pi*(campioni.d(conf.campione)/2)^2;

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

N   = length(result.indici_colpo.max_A); 
fs  = 52100;

% accelerazione/forza
lim_sup     = find(f > min(vertcat([PSD.fmax]))*1);
lim_inf     = find(f < 10);
f1          = f(lim_inf(end));

% definisco il vettore tempo
tempo       = 1:length(sng.F_filt(:,1));
tempo       = tempo/fs;
StopTime    = tempo(end);
simSng      = struct('V', zeros(size(sng.F)));
simPSD = struct('A',cell(1,N), 'A_fft',[]);

% massa
m3 = 13;
m2 = campioni.massa(conf.campione);
m1 = piastre.massa(conf.piastra);

% m3= m3 /3+2*m2 /3;
% m2= m2 /3+2*m1 /3;
% m1= m1/3;

M = [ m1+m2/3 m2 m3];

% rigidezza attesa
S_star=80*10^6;
 
% k attesa
k = S_star * s;

k1 = 40e9; %cost.elast. piastra carico grande [N/m]
k2 = k; % 2.8e+9; %10*k; %cost.elast.AC ref=4e6N/m
k3 = 21.6e+9; %cost.elast. base cemento [N/m];

% vettore delle K
K = [k1  k2  k3];

c1 = 0.00001*2*sqrt(K(1)*M(1)); 
c2 = 0.05*2*sqrt(K(2)*M(1));
c3 = 0.03*2*sqrt(K(3)*M(3)); 

% Vettore delle C
C = [c1 c2 c3];

% Simulo le velocità
[r, c]      = size(sng.F_filt);
for i=1:c
    %Definisco il timeseries Forza
    forza       = timeseries(sng.F_filt(:,i), tempo);
    
    % simulazione con simulink
    SimoutData = sim('Sim_sistema_fisico_ordine3.slx');
    
    % Acquisizione degli output di simulink
    simSng.V(:,i) = SimoutData.velocita.velocita.Data(1:end-1,1);
    simSng.A(:,i) = SimoutData.velocita.accelerazione.Data(1:end-1,1);
end

for i=1:N
    [simPSD(i).A, f] = periodogram(simSng(i).A, win_A, r, fs); %PSD Accelerazione [g^2]
    [simPSD(i).V, f] = periodogram(simSng(i).V, win_A, r, fs); %PSD Accelerazione [g^2]
end

% plot forza e velocità
figure

subplot(3,1,1)
plot(tempo, sng.F)

subplot(3,1,2)
plot(tempo, 20*log10(abs (simSng.V)))

subplot(3,1,3)
plot(tempo, 20*log10(abs (simSng.A)))

figure

subplot(2,1,1)
plot(f, 10*log10(abs (PSD.F)))
grid on, set(gca, 'XScale', 'log'), %xlim([ascissamin ascissamax])

subplot(2,1,2)
plot(f, 10*log10(abs (simPSD.A)))
grid on, set(gca, 'XScale', 'log'), %xlim([ascissamin ascissamax])

%%

f2 = f(lim_sup(1));
f_zoom = f1:0.2:f2;
prominanza = 0; %0.2;
distanza=20;
result.indici_colpo.max_A(1).locs = {1};
result.indici_colpo.max_A(1).A_norm_zoom = [];

for i = 1:N
    [~,cc] = size(PSD(i).F);
    [PSD(i).F_zoom, ~, pxxc] = periodogram(sng(i).F_filt, win_1, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Forza [N^2]
    [simPSD(i).A_zoom, ~, pyyc] = periodogram(simSng(i).A, win_A, f_zoom, fs, 'ConfidenceLevel',0.95); %PSD Accelerazione [g^2]
    
    figure(120+i), hold on
    for j = 1:cc
        % secondo minimo di K
        simPSD(i).K(:,j) = PSD(i).F(:,j) ./ simPSD(i).A(:,j) .* (2*pi*f).^4;
        plot (f, 10*log10(simPSD(i).K(:,j)), 'color', string(colore(i,:)));
        
        simPSD(i).K_zoom(:,j) = PSD(i).F_zoom(:,j) ./ simPSD(i).A_zoom(:,j) .* (1i*2*pi*f_zoom').^4;
        plot (f_zoom, 10*log10(simPSD(i).K_zoom(:,j)), 'color', string(colore(i+1,:)));
        
        [min_k, PosMin] = min(simPSD(i).K_zoom(lim_inf(end):lim_sup(1), j));
        fd = f_zoom(lim_inf(end) + PosMin );
        result.indici_colpo.min_K(i).fd(j) = fd;
        result.indici_colpo.min_K(i).K0(j) = (2*pi*fd)^2 * m;
        result.indici_colpo.min_K(i).S_star(j) = (2*pi*fd)^2 * m/s;
        result.indici_colpo.min_K(i).E(j)  = (2*pi*fd)^2 * m * h/s;
    end
    xl=xline(f(lim_inf(end)),'.',['Lim inf ',f(lim_inf(end)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    xl=xline(f(lim_sup(1)),'.',['Lim sup ',f(lim_sup(1)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    grid on
    set(gca, 'XScale', 'log')
    xlim([ascissamin ascissamax])

    % secondo massimo di A
    figure(130+i), hold on
    for j = 1:cc
        % seconfo massimo di A
        temp = simPSD(i).A(:,j) ./ PSD(i).F(:,j);
        plot (f, temp, 'color', string(colore(i,:)));
        
        temp_zoom = simPSD(i).A_zoom(:,j) ./ PSD(i).F_zoom(:,j);
        result.indici_colpo.max_A(i).A_norm_zoom(:,j) = temp_zoom;

        plot (f_zoom, temp_zoom, 'color', string(colore(i+1,:)));
        
        [pks, locs, wi, pr] = findpeaks(temp_zoom, f_zoom,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
        result.indici_colpo.max_A(i).locs{:,j} = locs;
        plot(result.indici_colpo.max_A(i).locs{j}, pks,'*')
        
        max_A = max(temp_zoom);
        if ~isempty(result.indici_colpo.max_A(i).locs{j})
            fd = result.indici_colpo.max_A(i).locs{j}(1);
        end        
        if length(result.indici_colpo.max_A(i).locs{j}) > 1
            f2=result.indici_colpo.max_A(i).locs{j}(2);
        end
        
        result.indici_colpo.max_A(i).fd(j) = fd;
        result.indici_colpo.max_A(i).f2(j) = f2;
        result.indici_colpo.max_A(i).K0(j) = (2*pi*fd)^2 * m;
        result.indici_colpo.max_A(i).S_star(j) = (2*pi*fd)^2 * m/s;
        result.indici_colpo.max_A(i).E(j)  = (2*pi*fd)^2 * m * h/s;
    end
    xl=xline(f(lim_inf(end)),'.',['Lim inf ',f(lim_inf(end)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    xl=xline(f(lim_sup(1)),'.',['Lim sup ',f(lim_sup(1)),' Hz']); xl.LabelVerticalAlignment = 'bottom';
    grid on
    set(gca, 'XScale', 'log')
    xlim([ascissamin f(lim_sup(1))]), %ylim([0 2*max_A])

end

figure (120), hold on,
for i = 1:N
    [~,cc] = size(PSD(i).F);
    for j = 1:cc
        %plot (max(sng(i).F_filt(:,j)), result.indici_colpo.min_K(i).S_star(j)/10^6, '*', 'color', string(colore(i,:)))
        plot (max(sng(i).F_filt(:,j)), result.indici_colpo.max_A(i).S_star(j)/10^6, 'o', 'color', string(colore(i+1,:)), 'lineWidth', 3)
    end
end
grid on, %ylim([140 200])
X = max(vertcat([ sng([1:N]).F_filt]));
Y = vertcat([result.indici_colpo.max_A([1:N]).S_star]) / 10^6;

modelfun = @(b,x)b(1)*exp(-(x-b(4))/b(2))+b(3)
% mdl = fitnlm(X, Y, modelfun, [100 50 10 0]); % Modello dei primi picchi
% [fit,R] = nlinfit(X, Y, modelfun, [200 150 150 0]);
% 
% X_exp = 1:max(X*1.5);
% plot (X_exp, feval(mdl,X_exp))
% title(['Rigidezza dinamica per unità di superficie vs picco della forza'])
% xlabel(['Force [N]'])
% ylabel(['Apparent Stiffness S [N/m^3]'])
% figure
% plot (X, mdl.WeightedResiduals.Studentized,'*')

k=1;
for i = 1:N
    [~,cc] = size(PSD(i).F);
    for j = 1 : cc
        if length(result.indici_colpo.max_A(i).locs{j}) > 1
            Y2(k) = (2*pi*result.indici_colpo.max_A(i).locs{j}(2))^2*m/s/10^6; %#ok<SAGROW>
            result.indici_colpo.max_A(i).locs{j}(2)
            X2(k) = max(sng(i).F_filt(:,j)); % forza
            k = k+1;
        end
    end
end

figure (120)
hold on
plot (X2, Y2, '*', 'color', string(colore(i,:)))
mdl2 = fitnlm(X2, Y2, modelfun, [100 50 10 0]); % Modello dei secondi picchi
[fit2,R2] = nlinfit(X2, Y2, modelfun, [200 150 150 0]);

X2_exp = 1:max(X2*1.5);
plot (X2_exp, feval(mdl2,X2_exp))
title(['Rigidezza dinamica per unità di superficie vs picco della forza'])
xlabel(['Force [N]'])
ylabel(['Apparent Stiffness S [N/m^3]'])
figure
plot (X2, mdl2.WeightedResiduals.Studentized,'*')

prominanza=1;
distanza=10;

[pks,locs,wi,p3] = findpeaks(temp_zoom, f_zoom,'MinPeakProminence',prominanza,'MinPeakDistance',distanza,'Annotate','extents');
    
% figure, hold on
% plot(f_zoom, temp_zoom)
% plot(locs, pks,'*')

%figure (201)

X3 = max(vertcat([sng.F_filt]));
[X3, Fsort] = sort(X3);
Y3 = (2*pi*f_zoom).^2*m/s;
% Y3 = f_zoom;
Z3 = vertcat([ result.indici_colpo.max_A.A_norm_zoom]);
Z3 = Z3(:,Fsort);
% mesh(X3, Y3, Z3)
% %ylim([0 6*10^8])
% % Create zlabel
% zlabel('Acc/Force','FontSize',18);
% 
% % Create ylabel
% ylabel('Corresponding Dynamic Stiffness [N/m^3]','FontSize',18);
% 
% % Create xlabel
% xlabel('Force [N]','FontSize',18);
%  zlim([0 10]);

fig1 = create_mash(X3, Y3/10^6, Z3);

ylim([0 200])
zlim([0 100])
xlim([0 300])
set(fig1.CurrentAxes ,'CLim',[7.14683052857818e-10 25],'Colormap',...
    [0 0 0.515625;0 0 0.533622382198953;0 0 0.551619764397906;0 0 0.569617146596859;0 0 0.587614528795812;0 0 0.605611910994764;0 0 0.623609293193717;0 0 0.64160667539267;0 0 0.659604057591623;0 0 0.677601439790576;0 0 0.695598821989529;0 0 0.713596204188482;0 0 0.731593586387435;0 0 0.749590968586388;0 0 0.76758835078534;0 0 0.785585732984293;0 0 0.803583115183246;0 0 0.821580497382199;0 0 0.839577879581152;0 0 0.857575261780105;0 0 0.875572643979058;0 0 0.89357002617801;0 0 0.911567408376963;0 0 0.929564790575916;0 0 0.947562172774869;0 0 0.965559554973822;0 0 0.983556937172775;0 0.00155431937172756 1;0 0.0195517015706804 1;0 0.0375490837696333 1;0 0.0555464659685861 1;0 0.073543848167539 1;0 0.0915412303664919 1;0 0.109538612565445 1;0 0.127535994764398 1;0 0.14553337696335 1;0 0.163530759162303 1;0 0.181528141361256 1;0 0.199525523560209 1;0 0.217522905759162 1;0 0.235520287958115 1;0 0.253517670157068 1;0 0.27151505235602 1;0 0.289512434554973 1;0 0.307509816753926 1;0 0.325507198952879 1;0 0.343504581151832 1;0 0.361501963350785 1;0 0.379499345549738 1;0 0.39749672774869 1;0 0.415494109947643 1;0 0.433491492146596 1;0 0.451488874345549 1;0 0.469486256544502 1;0 0.487483638743455 1;0 0.505481020942408 1;0 0.523478403141361 1;0 0.541475785340314 1;0 0.559473167539267 1;0 0.57747054973822 1;0 0.595467931937173 1;0 0.613465314136125 1;0 0.631462696335078 1;0 0.649460078534031 1;0 0.667457460732984 1;0 0.685454842931937 1;0 0.70345222513089 1;0 0.721449607329843 1;0 0.739446989528796 1;0 0.757444371727749 1;0 0.775441753926702 1;0 0.793439136125655 1;0 0.811436518324608 1;0 0.829433900523561 1;0 0.847431282722514 1;0 0.865428664921467 1;0 0.88342604712042 1;0 0.901423429319373 1;0 0.919420811518326 1;0 0.937418193717279 1;0 0.955415575916232 1;0 0.973412958115185 1;0 0.991410340314138 1;0.00940772251309085 1 0.990592277486909;0.0274051047120438 1 0.972594895287956;0.0454024869109968 1 0.954597513089003;0.0633998691099498 1 0.93660013089005;0.0813972513089027 1 0.918602748691097;0.0993946335078557 1 0.900605366492144;0.117392015706809 1 0.882607984293191;0.135389397905762 1 0.864610602094238;0.153386780104715 1 0.846613219895285;0.171384162303668 1 0.828615837696332;0.189381544502621 1 0.810618455497379;0.207378926701574 1 0.792621073298426;0.225376308900527 1 0.774623691099473;0.243373691099479 1 0.756626308900521;0.261371073298432 1 0.738628926701568;0.279368455497385 1 0.720631544502615;0.297365837696338 1 0.702634162303662;0.315363219895291 1 0.684636780104709;0.333360602094244 1 0.666639397905756;0.351357984293197 1 0.648642015706803;0.36935536649215 1 0.63064463350785;0.387352748691103 1 0.612647251308897;0.405350130890056 1 0.594649869109944;0.423347513089009 1 0.576652486910991;0.441344895287962 1 0.558655104712038;0.459342277486915 1 0.540657722513085;0.477339659685868 1 0.522660340314132;0.495337041884821 1 0.504662958115179;0.513334424083774 1 0.486665575916226;0.531331806282727 1 0.468668193717273;0.549329188481679 1 0.450670811518321;0.567326570680632 1 0.432673429319368;0.585323952879585 1 0.414676047120415;0.603321335078538 1 0.396678664921462;0.62131871727749 1 0.37868128272251;0.639316099476443 1 0.360683900523557;0.657313481675396 1 0.342686518324604;0.675310863874349 1 0.324689136125651;0.693308246073301 1 0.306691753926699;0.711305628272254 1 0.288694371727746;0.729303010471207 1 0.270696989528793;0.74730039267016 1 0.25269960732984;0.765297774869112 1 0.234702225130888;0.783295157068065 1 0.216704842931935;0.801292539267018 1 0.198707460732982;0.819289921465971 1 0.180710078534029;0.837287303664923 1 0.162712696335077;0.855284685863876 1 0.144715314136124;0.873282068062829 1 0.126717931937171;0.891279450261782 1 0.108720549738218;0.909276832460734 1 0.0907231675392657;0.927274214659687 1 0.0727257853403129;0.94527159685864 1 0.0547284031413602;0.963268979057593 1 0.0367310209424074;0.981266361256545 1 0.0187336387434547;0.999263743455498 1 0.000736256544501934;1 0.982738874345549 0;1 0.964741492146596 0;1 0.946744109947644 0;1 0.928746727748691 0;1 0.910749345549738 0;1 0.892751963350785 0;1 0.874754581151833 0;1 0.85675719895288 0;1 0.838759816753927 0;1 0.820762434554974 0;1 0.802765052356022 0;1 0.784767670157069 0;1 0.766770287958116 0;1 0.748772905759163 0;1 0.730775523560211 0;1 0.712778141361258 0;1 0.694780759162305 0;1 0.676783376963352 0;1 0.6587859947644 0;1 0.640788612565447 0;1 0.622791230366494 0;1 0.604793848167541 0;1 0.586796465968589 0;1 0.568799083769636 0;1 0.550801701570683 0;1 0.53280431937173 0;1 0.514806937172778 0;1 0.496809554973825 0;1 0.478812172774872 0;1 0.460814790575919 0;1 0.442817408376967 0;1 0.424820026178014 0;1 0.406822643979061 0;1 0.388825261780108 0;1 0.370827879581156 0;1 0.352830497382203 0;1 0.33483311518325 0;1 0.316835732984297 0;1 0.298838350785345 0;1 0.280840968586392 0;1 0.262843586387439 0;1 0.244846204188486 0;1 0.226848821989534 0;1 0.208851439790581 0;1 0.190854057591628 0;1 0.172856675392675 0;1 0.154859293193723 0;1 0.13686191099477 0;1 0.118864528795817 0;1 0.100867146596864 0;1 0.0828697643979117 0;1 0.0625 0;1 0.0538461538461537 0;1 0.0451923076923073 0;1 0.036538461538461 0;1 0.0278846153846146 0;1 0.0192307692307683 0;1 0.0105769230769219 0;1 0.00192307692307558 0;0.993269230769229 0 0;0.984615384615383 0 0;0.975961538461537 0 0;0.96730769230769 0 0;0.958653846153844 0 0;0.949999999999998 0 0;0.941346153846151 0 0;0.932692307692305 0 0;0.924038461538458 0 0;0.915384615384612 0 0;0.906730769230766 0 0;0.898076923076919 0 0;0.889423076923073 0 0;0.880769230769227 0 0;0.87211538461538 0 0;0.863461538461534 0 0;0.854807692307688 0 0;0.846153846153841 0 0;0.837499999999995 0 0;0.828846153846149 0 0;0.820192307692302 0 0;0.811538461538456 0 0;0.80288461538461 0 0;0.794230769230763 0 0;0.785576923076917 0 0;0.776923076923071 0 0;0.768269230769224 0 0;0.759615384615378 0 0;0.750961538461532 0 0;0.742307692307685 0 0;0.733653846153839 0 0;0.724999999999993 0 0;0.716346153846146 0 0;0.7076923076923 0 0;0.699038461538454 0 0;0.690384615384607 0 0;0.681730769230761 0 0;0.673076923076914 0 0;0.664423076923068 0 0;0.655769230769222 0 0;0.647115384615375 0 0;0.638461538461529 0 0;0.629807692307683 0 0;0.621153846153836 0 0;0.61249999999999 0 0;0.603846153846144 0 0;0.595192307692297 0 0;0.586538461538451 0 0;0.577884615384605 0 0;0.569230769230758 0 0;0.560576923076912 0 0;0.551923076923066 0 0;0.543269230769219 0 0;0.534615384615373 0 0;0.525961538461527 0 0;0.51730769230768 0 0;0.508653846153834 0 0;0.5 0 0]);

savefig('Waterfall_SvsF_sim.fig');

fig2 = create_mash(X3, Y3/10^6, 10*log10(Z3));
xlim([-0 300])
zlim([-10 20])
set(fig2.CurrentAxes ,'CLim',[-10 20],'Colormap',...
    [0 0 0.515625;0 0 0.533622382198953;0 0 0.551619764397906;0 0 0.569617146596859;0 0 0.587614528795812;0 0 0.605611910994764;0 0 0.623609293193717;0 0 0.64160667539267;0 0 0.659604057591623;0 0 0.677601439790576;0 0 0.695598821989529;0 0 0.713596204188482;0 0 0.731593586387435;0 0 0.749590968586388;0 0 0.76758835078534;0 0 0.785585732984293;0 0 0.803583115183246;0 0 0.821580497382199;0 0 0.839577879581152;0 0 0.857575261780105;0 0 0.875572643979058;0 0 0.89357002617801;0 0 0.911567408376963;0 0 0.929564790575916;0 0 0.947562172774869;0 0 0.965559554973822;0 0 0.983556937172775;0 0.00155431937172756 1;0 0.0195517015706804 1;0 0.0375490837696333 1;0 0.0555464659685861 1;0 0.073543848167539 1;0 0.0915412303664919 1;0 0.109538612565445 1;0 0.127535994764398 1;0 0.14553337696335 1;0 0.163530759162303 1;0 0.181528141361256 1;0 0.199525523560209 1;0 0.217522905759162 1;0 0.235520287958115 1;0 0.253517670157068 1;0 0.27151505235602 1;0 0.289512434554973 1;0 0.307509816753926 1;0 0.325507198952879 1;0 0.343504581151832 1;0 0.361501963350785 1;0 0.379499345549738 1;0 0.39749672774869 1;0 0.415494109947643 1;0 0.433491492146596 1;0 0.451488874345549 1;0 0.469486256544502 1;0 0.487483638743455 1;0 0.505481020942408 1;0 0.523478403141361 1;0 0.541475785340314 1;0 0.559473167539267 1;0 0.57747054973822 1;0 0.595467931937173 1;0 0.613465314136125 1;0 0.631462696335078 1;0 0.649460078534031 1;0 0.667457460732984 1;0 0.685454842931937 1;0 0.70345222513089 1;0 0.721449607329843 1;0 0.739446989528796 1;0 0.757444371727749 1;0 0.775441753926702 1;0 0.793439136125655 1;0 0.811436518324608 1;0 0.829433900523561 1;0 0.847431282722514 1;0 0.865428664921467 1;0 0.88342604712042 1;0 0.901423429319373 1;0 0.919420811518326 1;0 0.937418193717279 1;0 0.955415575916232 1;0 0.973412958115185 1;0 0.991410340314138 1;0.00940772251309085 1 0.990592277486909;0.0274051047120438 1 0.972594895287956;0.0454024869109968 1 0.954597513089003;0.0633998691099498 1 0.93660013089005;0.0813972513089027 1 0.918602748691097;0.0993946335078557 1 0.900605366492144;0.117392015706809 1 0.882607984293191;0.135389397905762 1 0.864610602094238;0.153386780104715 1 0.846613219895285;0.171384162303668 1 0.828615837696332;0.189381544502621 1 0.810618455497379;0.207378926701574 1 0.792621073298426;0.225376308900527 1 0.774623691099473;0.243373691099479 1 0.756626308900521;0.261371073298432 1 0.738628926701568;0.279368455497385 1 0.720631544502615;0.297365837696338 1 0.702634162303662;0.315363219895291 1 0.684636780104709;0.333360602094244 1 0.666639397905756;0.351357984293197 1 0.648642015706803;0.36935536649215 1 0.63064463350785;0.387352748691103 1 0.612647251308897;0.405350130890056 1 0.594649869109944;0.423347513089009 1 0.576652486910991;0.441344895287962 1 0.558655104712038;0.459342277486915 1 0.540657722513085;0.477339659685868 1 0.522660340314132;0.495337041884821 1 0.504662958115179;0.513334424083774 1 0.486665575916226;0.531331806282727 1 0.468668193717273;0.549329188481679 1 0.450670811518321;0.567326570680632 1 0.432673429319368;0.585323952879585 1 0.414676047120415;0.603321335078538 1 0.396678664921462;0.62131871727749 1 0.37868128272251;0.639316099476443 1 0.360683900523557;0.657313481675396 1 0.342686518324604;0.675310863874349 1 0.324689136125651;0.693308246073301 1 0.306691753926699;0.711305628272254 1 0.288694371727746;0.729303010471207 1 0.270696989528793;0.74730039267016 1 0.25269960732984;0.765297774869112 1 0.234702225130888;0.783295157068065 1 0.216704842931935;0.801292539267018 1 0.198707460732982;0.819289921465971 1 0.180710078534029;0.837287303664923 1 0.162712696335077;0.855284685863876 1 0.144715314136124;0.873282068062829 1 0.126717931937171;0.891279450261782 1 0.108720549738218;0.909276832460734 1 0.0907231675392657;0.927274214659687 1 0.0727257853403129;0.94527159685864 1 0.0547284031413602;0.963268979057593 1 0.0367310209424074;0.981266361256545 1 0.0187336387434547;0.999263743455498 1 0.000736256544501934;1 0.982738874345549 0;1 0.964741492146596 0;1 0.946744109947644 0;1 0.928746727748691 0;1 0.910749345549738 0;1 0.892751963350785 0;1 0.874754581151833 0;1 0.85675719895288 0;1 0.838759816753927 0;1 0.820762434554974 0;1 0.802765052356022 0;1 0.784767670157069 0;1 0.766770287958116 0;1 0.748772905759163 0;1 0.730775523560211 0;1 0.712778141361258 0;1 0.694780759162305 0;1 0.676783376963352 0;1 0.6587859947644 0;1 0.640788612565447 0;1 0.622791230366494 0;1 0.604793848167541 0;1 0.586796465968589 0;1 0.568799083769636 0;1 0.550801701570683 0;1 0.53280431937173 0;1 0.514806937172778 0;1 0.496809554973825 0;1 0.478812172774872 0;1 0.460814790575919 0;1 0.442817408376967 0;1 0.424820026178014 0;1 0.406822643979061 0;1 0.388825261780108 0;1 0.370827879581156 0;1 0.352830497382203 0;1 0.33483311518325 0;1 0.316835732984297 0;1 0.298838350785345 0;1 0.280840968586392 0;1 0.262843586387439 0;1 0.244846204188486 0;1 0.226848821989534 0;1 0.208851439790581 0;1 0.190854057591628 0;1 0.172856675392675 0;1 0.154859293193723 0;1 0.13686191099477 0;1 0.118864528795817 0;1 0.100867146596864 0;1 0.0828697643979117 0;1 0.0625 0;1 0.0538461538461537 0;1 0.0451923076923073 0;1 0.036538461538461 0;1 0.0278846153846146 0;1 0.0192307692307683 0;1 0.0105769230769219 0;1 0.00192307692307558 0;0.993269230769229 0 0;0.984615384615383 0 0;0.975961538461537 0 0;0.96730769230769 0 0;0.958653846153844 0 0;0.949999999999998 0 0;0.941346153846151 0 0;0.932692307692305 0 0;0.924038461538458 0 0;0.915384615384612 0 0;0.906730769230766 0 0;0.898076923076919 0 0;0.889423076923073 0 0;0.880769230769227 0 0;0.87211538461538 0 0;0.863461538461534 0 0;0.854807692307688 0 0;0.846153846153841 0 0;0.837499999999995 0 0;0.828846153846149 0 0;0.820192307692302 0 0;0.811538461538456 0 0;0.80288461538461 0 0;0.794230769230763 0 0;0.785576923076917 0 0;0.776923076923071 0 0;0.768269230769224 0 0;0.759615384615378 0 0;0.750961538461532 0 0;0.742307692307685 0 0;0.733653846153839 0 0;0.724999999999993 0 0;0.716346153846146 0 0;0.7076923076923 0 0;0.699038461538454 0 0;0.690384615384607 0 0;0.681730769230761 0 0;0.673076923076914 0 0;0.664423076923068 0 0;0.655769230769222 0 0;0.647115384615375 0 0;0.638461538461529 0 0;0.629807692307683 0 0;0.621153846153836 0 0;0.61249999999999 0 0;0.603846153846144 0 0;0.595192307692297 0 0;0.586538461538451 0 0;0.577884615384605 0 0;0.569230769230758 0 0;0.560576923076912 0 0;0.551923076923066 0 0;0.543269230769219 0 0;0.534615384615373 0 0;0.525961538461527 0 0;0.51730769230768 0 0;0.508653846153834 0 0;0.5 0 0]);




