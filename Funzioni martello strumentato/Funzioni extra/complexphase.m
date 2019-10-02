function [val] = complexphase(num)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
val = atan(imag(num)./real(num));
end

