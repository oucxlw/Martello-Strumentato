%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Modello matematico per i campioni di piccole dimensioni
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

clear all
%close all

    f=[10:80000];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           f=2:800000;
    f=2*pi.*f;


for a=2:0.5:4

    m1=0.6;
    m2=1.5;
%     s=0.10*0.097;

    k1=0.98*1e9*0.04;
    k2=300e9*0.008;

    c1=a*1;
    c2=0.03;

    k12=-k2;
    k11=k1+k2;
    k22=k2;
    k21=k12;

    A = (-k21+ 1i*f*c2) ./ (-m1*f.^2 + k11 - 1i*f*(c1+c2));
    A2 = (-k21- 1i*f*c2) ./ (m1*f.^2 + k11 - 1i*f*(c1+c2));
    displacement = -f.^2*m2 + k22*A- 1i.*f.*c2.*(A-1);
    Dmass   = (-f.^2*m2 + k22+k12*A2 - 1i.*f.*c2.*(A-1)) ./  f.^2;
    Dmass2 = (+f.^2*m2 + k22+ k12*A2 +1i.*f.*c2.*(A2-1)) ./  f.^2;

    figure (202)
    set(gca, 'XScale', 'log')
    hold on
    plot (f, abs(Dmass))
    %plot (f, abs(Dmass2))

    xlim([20, 8000])
end
set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
hold off

%%
puppa