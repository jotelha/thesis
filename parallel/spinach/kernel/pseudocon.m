% Calculates the pseudocontact shift field of a point metal centre 
% and adds the corresponding shifts to the chemical shielding tensors
% of all nuclei. Input parameters:
%
%     metal_coords  -  [X Y Z] coordinates of the metal, in Angstroms
%     chi           -  magnetic susceptibility tensor of the metal,
%                      as a 3x3 matrix in the molecular reference,
%                      frame, in units of cubic Angstroms.
%
% ilya.kuprov@oerc.ox.ac.uk
% gareth.charnock@oerc.ox.ac.uk

function spin_system=pseudocon(spin_system,metal_coords,chi)

% Isolate the second-rank component of chi
[~,~,rank2]=mat2sphten(chi); chi=sphten2mat([],[],rank2);

% Loop over the spins in the system
for n=1:spin_system.comp.nspins
    
    % Bomb out if the system contains electrons
    if strcmp(spin_system.comp.isotopes{n}(1),'E')
        error('pseudocon: PCS calculations are only available for NMR systems.');
    end
    
    % Compute the PCS for the current nucleus
    x=(spin_system.inter.coordinates(n,1)-metal_coords(1));
    y=(spin_system.inter.coordinates(n,2)-metal_coords(2));
    z=(spin_system.inter.coordinates(n,3)-metal_coords(3));
    sigma=(1/(12*pi))*((2*x^2-y^2-z^2)*chi(1,1)+...
                       (2*y^2-x^2-z^2)*chi(2,2)+...
                       (2*z^2-x^2-y^2)*chi(3,3)+...
                        6*x*y*chi(1,2)+...
                        6*x*z*chi(1,3)+...
                        6*y*z*chi(2,3))./((x^2+y^2+z^2).^(5/2));
                      
    % Add the PCS contribution to the Zeeman tensor
    spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+...
                                       sigma*spin_system.inter.basefrq{n}*eye(3);
                                   
end

end

% It's not true I had nothing on. I had the radio on.
%
% Marilyn Monroe

