% Isotropic singlet-singlet MARY experiment with exponential recombination
% function.
%
% fields:   row vector of field values in Tesla
% kinetics: row vector of singlet recombination rate constants in Hz
% L_hfc:    optional Liouvillian, generated from scratch if omitted
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function M=mary(spin_system,fields,kinetics,L_hfc)

% Show the banner
banner(spin_system,'sequence_banner');

% Get the fields into angular frequencies
gamma=spin('E'); fields=fields*gamma;

% Get the two-electron singlet state
electrons=find(strcmp(spin_system.comp.isotopes,'E'));
if length(electrons)~=2
    error('mary: number of electrons not equal to two.');
end
S=singlet(spin_system,electrons(1),electrons(2));

% Generate the field operators
report(spin_system,'mary: generating the field operators...');
E1z=operator(spin_system,'Lz',electrons(1));
E2z=operator(spin_system,'Lz',electrons(2));

% Get the Liouvillian
switch nargin
    case 3
        report(spin_system,'mary: building the Liouvillian...');
        L_hfc=h_superop(secularity(spin_system,'lowfield'))+1i*r_superop(spin_system)+1i*k_superop(spin_system);
    case 4
        report(spin_system,'mary: using the Liouvillian as supplied...');
end

% Run the state space restriction
report(spin_system,'mary: trying to reduce the problem dimension...');
projectors=reduce(spin_system,L_hfc+mean(abs(fields))*(E1z+E2z),S);

% Preallocate the result array
M=zeros(length(kinetics),length(fields));

% Loop over the independent subspaces
for subs=1:length(projectors)
    
    % Project into the current subspace
    L_hfc_subs=projectors{subs}'*L_hfc*projectors{subs};
    E1z_subs=projectors{subs}'*E1z*projectors{subs};
    E2z_subs=projectors{subs}'*E2z*projectors{subs};
    S_subs=projectors{subs}'*S;
    
    % Loop over the rate constants and fields
    for n=1:length(kinetics)
        for m=1:length(fields)
            
            % Assemble the Liouvillian
            L=L_hfc_subs+fields(m)*(E1z_subs+E2z_subs);
            
            % Compute MARY
            M(n,m)=M(n,m)+kinetics(n)*real(S_subs'*((1i*L+kinetics(n)*speye(size(L)))\S_subs));
            
        end
    end
    
end

end