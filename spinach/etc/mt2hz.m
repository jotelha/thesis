% Converts hyperfine couplings from milliTesla to Hz (linear frequency).
%
% ilya.kuprov@oerc.ox.ac.uk

function hfc_hz=mt2hz(hfc_mt)

    hfc_hz=1e7*2.802495365*hfc_mt;
    
end

% You know that I write slowly. This is chiefly because I am never
% satisfied until I have said as much as possible in a few words, and
% writing briefly takes far more time than writing at length.
%
% Carl Friedrich Gauss 