%%
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Modello matematico per i campioni di piccole dimensioni
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

clear all
%close all

    f=logspace(0 ,4,10000)%[10:80000];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           f=2:800000;
    f=2*pi.*f;


for a=0.0001

    m1=30;
    m2=1.81;
%   s=0.10*0.097;



    k1=1;
    k2=1.5e6;

    c1=1;
    c2=1;

    c2=c2*sqrt(k2*m2)
    
    c1=c1*sqrt(k1*(m1+m2))
    
    k12=-k2;
    k11=k1+k2;
    k22=k2;
    k21=k12;

    A = (-k21+ 1i*f*c2) ./ (-m1*f.^2 + k11 - 1i*f*(c1+c2));
    A2 = (-k21- 1i*f*c2) ./ (m1*f.^2 + k11 - 1i*f*(c1+c2));
    displacement = -f.^2*m2 + k22*A- 1i.*f.*c2.*(A-1);
    Dmass   = (-f.^2*m2 + k22+k12*A2 - 1i.*f.*c2.*(A-1)) ;
    Dmass2 = (+f.^2*m2 + k22+ k12*A2 +1i.*f.*c2.*(A2-1)) ;

%     figure (202)
    set(gca, 'XScale', 'log')
    hold on
    plot (f/(2*pi), 20*log10(abs(Dmass)))
    %plot (f, abs(Dmass2))

    xlim([20, 1e4])
end
set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
hold off

%%
puppa