% A simulation of solid effect DNP in a system with five
% nuclei using third-order Krylov-Bogolyubov averaging.
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

function dnp_solid_kb_third_order()

% Spin system
sys.magnet=3.4;
sys.isotopes={'1H','1H','1H','1H','1H','E'};

% Basis: all coherences of order 0, +1 and -1
bas.mode='complete';
bas.projections=[-1 0 1];

% Relaxation theory parameters
inter.equilibrium='thermal';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Nuclear Zeeman frequency (rad/s)
omega=2*pi*144.7e6;

% Microwave irradiation power (rad/s)
M=2*pi*1e6;

% Time stepping
time_step=0.3;
n_steps=1e3;

% Empirical relaxation rates (Hz)
r1e=1e3; r2e=1e6;
r1n=1e-2; r2n=1e0;

% Parameters of hyperfine interactions (random)
A(1)=-2*pi*50851.03; A(2)=-2*pi*56038.01;
A(3)=-2*pi*28037.63; A(4)=-2*pi*16603.86;
A(5)=2*pi*32829.34;
B(1)=-2*pi*11066.93; B(2)=0;
B(3)=-2*pi*6101.948; B(4)=0;
B(5)=2*pi*3531.22;

% Parameters of nuclear dipolar interactions (random)
d(1,2)=2*pi*29; d(1,3)=2*pi*7; d(1,4)=2*pi*14;
d(1,5)=-2*pi*18; d(2,3)=2*pi*20; d(2,4)=2*pi*5;
d(2,5)=-2*pi*15; d(3,4)=2*pi*11; d(3,5)=-2*pi*6;
d(4,5)=-2*pi*4;

tic,

% Locate the electron and the nuclei
electron_position=find(strcmp(spin_system.comp.isotopes,'E'));
nuclear_positions=find(~strcmp(spin_system.comp.isotopes,'E'));

% Get electron spin operators
Sp=operator(spin_system,'L+',electron_position);
Sz=operator(spin_system,'Lz',electron_position);

% Preallocate Hamiltonian components
matrix_dim=spin_system.bas.nstates;
V0=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
VpSm=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
VpSz=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
V=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
HD=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
RI=spalloc(matrix_dim,matrix_dim,10*matrix_dim);

% Loop over the nuclear spins
for n=1:numel(nuclear_positions)
    
    % Get electron-nuclear two-spin operators
    IzSz=operator(spin_system,{'Lz','Lz'},{nuclear_positions(n),electron_position});
    IpSm=operator(spin_system,{'L+','L-'},{nuclear_positions(n),electron_position});
    IpSz=operator(spin_system,{'L+','Lz'},{nuclear_positions(n),electron_position});
    
    % Get nuclear spin operators
    Iz=operator(spin_system,'Lz',nuclear_positions(n));
    Ip=operator(spin_system,'L+',nuclear_positions(n));
    
    % Update the coherent superoperators (see the paper)
    V0=V0+A(n)*IzSz;
    VpSm=VpSm+0.5*B(n)*IpSm;
    VpSz=VpSz+0.5*B(n)*IpSz;
    V=V+0.5*(abs(B(n)))^2*Iz;
    
    % Update the relaxation superoperator (see the paper)
    RI=RI+0.25*r1n*(Ip'*Ip+Ip*Ip')+(r2n-0.5*r1n)*Iz*Iz;
    
end

% Loop over pairs of nuclear spins
for n=1:numel(nuclear_positions)
    for m=(n+1):numel(nuclear_positions)
        
        % Get nuclear two-spin operators
        IzIz=operator(spin_system,{'Lz','Lz'},{nuclear_positions(n),nuclear_positions(m)});
        IpIm=operator(spin_system,{'L+','L-'},{nuclear_positions(n),nuclear_positions(m)});
                
        % Update the dipolar Hamiltonian superoperator
        HD=HD+d(n,m)*(2*IzIz-0.5*(IpIm'+IpIm));
        
    end
end

% Compute the relaxation superoperator
R=-0.25*r1e*(Sp'*Sp+Sp*Sp')-(r2e-0.5*r1e)*Sz*Sz-RI;

% Run third-order Krylov-Bogolyubov averaging (see the paper)
lp=0.5*M*Sp+VpSz;
l0=-V0-1i*R;
clml0=lp'*l0-l0*lp';
clpl0=lp*l0-l0*lp;
H1=HD+V0;
H2=0.25*(V+2*M^2*Sz-2*M*(VpSm+VpSm'))/omega;
H3=(lp*clml0+lp'*clpl0)/omega^2;
H=H1+H2+H3;

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
       state(spin_system,'Lz',4)...
       state(spin_system,'Lz',5)...
       state(spin_system,'Lz',6)];

% Simulation
answer=evolution(spin_system,L,coils,rho,time_step,n_steps,'multichannel');

% Plotting
plot(real(answer'));

end

