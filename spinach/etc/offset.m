% Generates the offset superoperator.
%
% When added to the main Liouvillian, the offset superoperator shifts the
% master frame of a given set of spins by a given amount (in Hz). Syntax:
%
%                   A=offset(spin_system,spins,offset_value)
%
% where 'spins' may be set to 'all', '1H', '13C' etc.
%
% Ilya Kuprov (ilya.kuprov@oerc.ox.ac.uk)
% Konstantin Pervushin (p@trosy.com)

function A=offset(spin_system,spins,offset_value)

    A=2*pi*offset_value*operator(spin_system,'Lz',spins);

end

% Physics is becoming too difficult for the physicists.
% 
% David Hilbert 