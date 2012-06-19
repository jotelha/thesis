% A simulation of solid effect DNP in a system with three
% nuclei using fourth-order Krylov-Bogolyubov averaging.
%
% The electron-nucleus zero-quantum transition is irradi-
% ated at omega_e+omega_n.
%
% Empirical relaxation model is used, as described in our 
% PCCP paper at http://xlink.rsc.org/?doi=C2CP23233B
%
% Alexander Karabanov (Nottingham)
% Anniek van der Drift (Nottingham)
% Walter Kockenberger (Nottingham)
% Luke Edwards (Oxford)
% Ilya Kuprov (Oxford)

function dnp_solid_kb_fourth_order()

% Spin system
sys.magnet=3.4;
sys.isotopes={'1H', '1H', '1H', 'E'};

% Complete basis set
bas.mode='complete';

% Relaxation theory parameters
inter.equilibrium='thermal';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Nuclear Zeeman frequency (rad/s)
omega=2*pi*144.7e6;

% Microwave irradiation power (rad/s)
M=2*pi*2.5e6;

% Time stepping
time_step=1.5e-3;
n_steps=1e3;

% Empirical relaxation rates (Hz)
r1e=1e3; r2e=1e6;
r1n=1e-2; r2n=1e0;

% Locate the electron and the nuclei
electron_position=find(strcmp(spin_system.comp.isotopes,'E'));
nuclear_positions=find(~strcmp(spin_system.comp.isotopes,'E'));

% Get electron spin operators
Sp=operator(spin_system,'L+',electron_position);
Sm=operator(spin_system,'L-',electron_position);
Sz=operator(spin_system,'Lz',electron_position);

% Preallocate Hamiltonian components
matrix_dim=spin_system.bas.nstates;
V0=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
VpSm=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
VmSp=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
VpSz=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
VmSz=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
V=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
HD=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
RI=spalloc(matrix_dim,matrix_dim,10*matrix_dim);

% Loop over the nuclear spins
for n=1:spin_system.comp.nspins-1
    
    % Get electron-nuclear two-spin operators
    IzSz=operator(spin_system,{'Lz','Lz'},{nuclear_positions(n),electron_position});
    ImSp=operator(spin_system,{'L-','L+'},{nuclear_positions(n),electron_position});
    IpSm=operator(spin_system,{'L+','L-'},{nuclear_positions(n),electron_position});
    IpSz=operator(spin_system,{'L+','Lz'},{nuclear_positions(n),electron_position});
    ImSz=operator(spin_system,{'L-','Lz'},{nuclear_positions(n),electron_position});
    
    % Get nuclear spin operators
    Iz=operator(spin_system,'Lz',nuclear_positions(n));
    Im=operator(spin_system,'L-',nuclear_positions(n));
    Ip=operator(spin_system,'L+',nuclear_positions(n));
    
    % Set the interaction amplitudes (random)
    A=1e6/n*2*pi; B=A/2;
    
    % Update the coherent superoperators (see the paper)
    V0=V0+A*IzSz;
    VpSm=VpSm+0.5*B*IpSm;
    VmSp=VmSp+0.5*B*ImSp;
    VpSz=VpSz+0.5*B*IpSz;
    VmSz=VmSz+0.5*B*ImSz;
    V=V+0.5*B^2*Iz;
    
    % Update the relaxation superoperator (see the paper)
    RI=RI+0.25*r1n*(Im*Ip+Ip*Im)+(r2n-0.5*r1n)*Iz*Iz;
    
end

% Loop over pairs of nuclear spins
for n=1:numel(nuclear_positions)
    for m=(n+1):numel(nuclear_positions)
        
        % Set the dipolar interaction parameter (random)
        d=1e3*2*pi/((n-m)^2);
        
        % Get nuclear two-spin operators
        IzIz=operator(spin_system,{'Lz','Lz'},{nuclear_positions(n),nuclear_positions(m)});
        IpIm=operator(spin_system,{'L+','L-'},{nuclear_positions(n),nuclear_positions(m)});
        ImIp=operator(spin_system,{'L-','L+'},{nuclear_positions(n),nuclear_positions(m)});
        
        % Update the dipolar Hamiltonian superoperator
        HD=HD+d*(2*IzIz-0.5*(ImIp+IpIm));
        
    end
end

% Compute the relaxation superoperator
R=-0.25*r1e*(Sm*Sp+Sp*Sm)-(r2e-0.5*r1e)*Sz*Sz-RI;

% Run fourth-order Krylov-Bogolyubov averaging (see the paper)
lp=0.5*M*Sp+VpSz; lm=0.5*M*Sm+VmSz; l0=HD+V0-1i*R;
H=HD+V0+0.25*(V+2*M^2*Sz-2*M*(VpSm+VmSp))/omega+...
             (lp*(lm*l0-l0*lm)+lm*(lp*l0-l0*lp))/omega^2-...
             (0.5*(lm*lm*lp*lp-lp*lp*lm*lm)+(lp*lm+lm*lp)*(lp*lm-lm*lp)+...
              lp*(l0*(lm*l0-l0*lm)-(lm*l0-l0*lm)*l0)-...
              lm*(l0*(lp*l0-l0*lp)-(lp*l0-l0*lp)*l0))/omega^3;

% Set the thermal equilibrium state to electron Lz
sigma_th=state(spin_system,'Lz',electron_position);
R(:,1)=-R*sigma_th;

% Assemble the Liouvillian superoperator
L=H+1i*R;

% Set the initial state to thermal equilibrium
rho=sigma_th;

% Set the detection states to Lz on every spin
coils=[state(spin_system,'Lz',1)...
       state(spin_system,'Lz',2)...
       state(spin_system,'Lz',3)...
       state(spin_system,'Lz',4)];

% Simulation
answer=evolution(spin_system,L,coils,rho,time_step,n_steps,'multichannel');

% Plotting
plot(real(answer'));

end

