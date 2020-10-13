
function [dat_cut] = taglio_completo(dat, inizio_misura_1, fine_misura_1, inizio_misura_2, fine_misura_2)
for i = 1:3
    dat_cut(i).transiti = taglio(dat(i).transiti, inizio_misura_1, fine_misura_1);
    dat_cut(i).livelli = taglio(dat(i).livelli , inizio_misura_1, fine_misura_1);
    
    dat_cut(i).durations_transiti = seconds(dat_cut(i).transiti.Time(end) - dat_cut(i).transiti.Time(1));
    dat_cut(i).durations_livelli = seconds(dat_cut(i).livelli.Time(end) - dat_cut(i).livelli.Time(1));

end

for i = 4:6
    dat_cut(i).transiti = taglio(dat(i).transiti, inizio_misura_2, fine_misura_2);
    dat_cut(i).livelli = taglio(dat(i).livelli , inizio_misura_2, fine_misura_2);
    
    dat_cut(i).durations_transiti = seconds(dat_cut(i).transiti.Time(end) - dat_cut(i).transiti.Time(1));
    dat_cut(i).durations_livelli = seconds(dat_cut(i).livelli.Time(end) - dat_cut(i).livelli.Time(1));

end
dat_cut(7).transiti = taglio(dat(i).transiti, inizio_misura_2, fine_misura_2);
dat_cut(8).transiti = taglio(dat(i).transiti, inizio_misura_2, fine_misura_2);

dat_cut(7).durations_transiti = seconds(dat_cut(7).transiti.Time(end) - dat_cut(7).transiti.Time(1));
dat_cut(8).durations_transiti = seconds(dat_cut(8).transiti.Time(end) - dat_cut(8).transiti.Time(1));

end