load dati.mat
data(:,1)=data(:,1)*2000;

data(:,2)=data(:,2)*500;
[pks,picchi,w,p] = findpeaks(data(:,1),fs,'MinPeakProminence',40,'MinPeakDistance',0.5,'Annotate','extents');

length(picchi);

t=(0:length(data)-1)/fs;
figure, hold on
plot(t,data(:,1))
plot(picchi,pks,'*')