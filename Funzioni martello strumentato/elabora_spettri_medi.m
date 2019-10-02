function [PSD_Fav, PSD_Aav,FFT_Fav, FFT_Aav,Cxy,f]=elabora_spettri_medi(Forza,Accelerazione, soglia, delay, inizio, fine, L_pre, L_coda, filt_doppi,window_F,window_A,bandwidth,fs)


[picchi,n_picchi] = trovacolpi(Forza, soglia, delay, inizio, fine);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Definizione delle matrici (selezione dei segnali)
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[F, pos] = creamatriceforza (Forza, picchi,n_picchi, L_pre, L_coda, filt_doppi, fs);
[A] = creamatriceaccelerazione (Accelerazione, pos, L_pre, L_coda, fs);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Finestratura e Normalizzazione
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[F_filt] = finestra_forza (F, window_F, fs);
[A_filt] = finestra_accelerazione (A, window_A, fs);

%<<<<<<<<<<<<<<<<<<<<
% Calcolo delle FRFs
%<<<<<<<<<<<<<<<<<<<<
L = length(F_filt(:,1));
[PSD_F, f]= periodogram(F_filt, [], L, fs); %PSD Forza [N^2]
[PSD_A, f]= periodogram(A_filt, [], L, fs); %PSD Accelerazione [g^2]

%<<<<<<<<<<<<<<<<<<<<
% Filtraggio Banda
%<<<<<<<<<<<<<<<<<<<<
[r,c]=size(PSD_F);
tagli=[];
scarti=0;

for jj=1:(c-scarti)
f0=find(f>10,1);
  fmax=find(PSD_F(f0:end, jj-scarti)<((PSD_F(f0, jj-scarti)/10)),1); 
     fmax=f(fmax+f0);
     if  fmax<bandwidth
        PSD_F(:,jj-scarti)=[];
        tagli=[tagli; jj];
        scarti=scarti+1;
     end
end

PSD_A(:,tagli)=[];

% Matrici F e An filtrate con il bandwidth
F_filt2 = F_filt ;
A_filt2 = A_filt;
F_filt2(:,tagli)=[];
A_filt2(:,tagli)=[];

F_filtall = reshape(F_filt2, [],1);
A_filtall = reshape(A_filt2, [],1);
[r,c]=size(F_filt2);
[Cxy,f] = mscohere(F_filtall, A_filtall, round(length(F_filtall)./c),[],L,fs);

%Calcolo e memorizzo Spettri in fft di F, A, V e D


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% calcolo la media degli spettri
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%PSD
PSD_Fav = mean(sqrt(PSD_F), 2);
PSD_Aav = mean(sqrt(PSD_A), 2);
FFT_Fav = mean( fft (F_filt2,  [], 1), 2); %FFT forza
FFT_Aav = mean( fft (A_filt2, [], 1), 2); %FFT accelerazione
end