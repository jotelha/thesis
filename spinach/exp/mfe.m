% Magnetic field effect on the single state dynamics.
%
% field:    magnetic field in Tesla
% L_hfc:    optional Liouvillian, generated from scratch if omitted
%
% ilya.kuprov@oerc.ox.ac.uk

function fid=mfe(spin_system,field,timestep,nsteps)

% Show the banner
banner(spin_system,'sequence_banner');

% Get the fields into angular frequencies
gamma=spin('E'); field=field*gamma;

% Get the two-electron singlet state
if numel(spin_system.chem.rp_spins)
    S=singlet(spin_system,spin_system.chem.rp_spins(1),spin_system.chem.rp_spins(2));
else
    error('mfe: radical pair kinetics model should be specified in inter.chem.');
end

% Generate the field operators
E1z=operator(spin_system,'Lz',spin_system.chem.rp_spins(1));
E2z=operator(spin_system,'Lz',spin_system.chem.rp_spins(2));

% Get the Liouvillian
L=field*(E1z+E2z)+h_superop(secularity(spin_system,'lowfield'))+1i*r_superop(spin_system)+1i*k_superop(spin_system);

% Run the simulation
fid=evolution(spin_system,L,S,S,timestep,nsteps,'observable');

end