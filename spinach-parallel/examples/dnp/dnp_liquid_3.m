% Overhauser type DNP in liquid phase at room temperature, using a
% single inversion pulse on the electron ESR signal.
%
% Ilya Kuprov (Oxford)
% Luke Edwards (Oxford)
% Alexander Karabanov (Nottingham)
% Anniek van der Drift (Nottingham)
% Walter Kockenberger (Nottingham)

function dnp_liquid_3()

% Spin system
sys.magnet=3.4;
sys.isotopes={'1H','1H','E'};

% Zeeman interactions
inter.zeeman.matrix={[5 0 0; 0 5 0; 0 0 5]
                     [5 0 0; 0 5 0; 0 0 5]
                     [2.0023 0 0; 0 2.0025 0; 0 0 2.0027]};

% Coordinates
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 2.0 0.0]
                   [0.0 0.0 1.5]};
               
% Complete basis set
bas.mode='complete';

% Relaxation theory
inter.relaxation='redfield';
inter.equilibrium='thermal';
inter.temperature=298;
inter.tau_c=10e-12;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial state - thermal equilibrium
rho=equilibrium(spin_system);

% Detection state - Lz on the protons
coil_a=statevec(spin_system,[2 0 0]);
coil_b=statevec(spin_system,[0 2 0]);

% Static Liouvillian superoperator
L=h_superop(secularity(spin_system,'nmr'));

% Electron control operator
Lp=operator(spin_system,'L+','E');

% Relaxation superoperator
R=r_superop(spin_system);

% Pi pulse on electron
rho=step(spin_system,(Lp+Lp')/2,rho,pi);

% Evolution (10000 steps, 0.1 microseconds each)
answer_a=evolution(spin_system,L+1i*R,coil_a,rho,1e-7,10000,'observable');
answer_b=evolution(spin_system,L+1i*R,coil_b,rho,1e-7,10000,'observable');

% Plotting
hold on;
plot(real(answer_a'));
plot(real(answer_b'));

end

