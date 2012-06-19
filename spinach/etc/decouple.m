% Basic analytical decoupler. 
%
% Disconnects a user-specified spin subsystem from the rest of the system
% by zeroing the corresponding rows and columns of the Liouvillian. The
% disconnected subsystem would stay frozen in time. Syntax:
%
%   [L,rho]=decouple(spin_system,L,rho,spins)
%
%     spins: spins to be decoupled, specified either by name ('13C', etc.)
%            or by number ([1 2 3], etc.).
%     L:     Liouvillian superoperator
%     rho:   state vector or a horizontal stack thereof
%
% ilya.kuprov@oerc.ox.ac.uk

function [L,rho]=decouple(spin_system,L,rho,spins)

% Find the nuclei to be decoupled
if isnumeric(spins)
    dec_mask=false(1,spin_system.comp.nspins); dec_mask(spins)=true;
else
    dec_mask=strcmp(spin_system.comp.isotopes,spins);
end

% Get the list of states to be disconnected
zero_mask=(sum(spin_system.bas.basis(:,dec_mask),2)~=0);

% Zero the corresponding rows and columns of the Liouvillian
L(zero_mask,:)=0; L(:,zero_mask')=0;

% Zero the corresponding rows of the state vector stack
if nargout==2, rho(zero_mask,:)=0; end

end

% Blessed are the hearts that can bend; they shall never be broken.
%
% Albert Camus 