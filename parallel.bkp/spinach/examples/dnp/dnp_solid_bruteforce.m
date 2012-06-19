% A brute-force simulation of solid effect DNP in a system 
% with three nuclei. The electron-nucleus zero-quantum tran-
% sition is irradiated at omega_e+omega_n.
%
% Empirical relaxation model is used, as described in our 
% PCCP paper at http://xlink.rsc.org/?doi=C2CP23233B
%
% Alexander Karabanov (Nottingham)
% Anniek van der Drift (Nottingham)
% Walter Kockenberger (Nottingham)
% Luke Edwards (Oxford)
% Ilya Kuprov (Oxford)

function dnp_solid_bruteforce()

% Spin system
sys.magnet=3.4;
sys.isotopes={'E','1H','1H','1H'};

% Coordinates
sys.tols.prox_cutoff=Inf;
inter.coordinates={[0 0 0]; [0 0 7]; [0 0 10]; [0 0 14]};

% Relaxation theory
inter.relaxation='t1_t2';
inter.r1_rates=[1e3 10 1.0 0.1];
inter.r2_rates=[1e6 100 10 1.0];
inter.temperature=4.2;
inter.equilibrium='thermal';

% Microwave irradiation power and offset
mw_power=2*pi*1.5e6;
mw_offset=2*pi*144.7e6;

% Time stepping
time_step=0.01;
n_steps=300;

% Basis: all coherences of order 0, +1 and -1
bas.mode='complete';
bas.projections=[-1 0 1];

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Electron rotating frame Hamiltonian superoperator
[h_iso,rotational_basis]=h_superop(secularity(spin_system,'solid_effect'));
H=h_iso+orientation(rotational_basis,[0 pi/4 0]);

% Offset term
Ez=operator(spin_system,'Lz','E');
H=H+mw_offset*Ez;

% Microwave irradiation term
Ep=operator(spin_system,'L+','E');
H=H+mw_power*(Ep+Ep')/2;

% Relaxation superoperator
L=H+1i*r_superop(spin_system);

% Simulation
rho=equilibrium(spin_system);
coils=[state(spin_system,'Lz',1)...
       state(spin_system,'Lz',2)...
       state(spin_system,'Lz',3)...
       state(spin_system,'Lz',4)];
answer=evolution(spin_system,L,coils,rho,time_step,n_steps,'multichannel');

% Plotting
plot(real(answer'));

end

