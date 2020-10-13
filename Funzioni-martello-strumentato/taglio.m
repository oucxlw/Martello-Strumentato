function [dati] = taglio(dati, in1, out1)

cut_in = find(dati.Time > in1, 1); %&& dati.Time.milliseconds==in1.milliseconds);
cut_out = find(dati.Time < out1); %&& dati.Time.milliseconds==out1.milliseconds);
cut_out = cut_out(end);

if dati.Time(end) ~= out1 %|| dati.Time(end).milliseconds ~= out1.milliseconds
    dati((cut_out+1:end), :) = [];
end

if dati.Time(1) ~= in1
    dati((1:cut_in-1), :) = [];
end

end