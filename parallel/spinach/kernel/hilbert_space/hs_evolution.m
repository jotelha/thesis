% Time evolution function in Hilbert space. Arguments:
%
%       H           - Hamiltonian matrix, must be Hermitian
%       coil        - observable operator (if any)
%       rho         - initial density matrix
%       timestep    - duration of a single time step (seconds)
%       nsteps      - number of steps to take
%       output      - output type: 'final' computes the final density matrix using
%                     the split algorithm (called B1 in the JMR paper), 'final_dt'
%                     does the same using the double transpose algorithm (called B2
%                     in the JMR paper), 'observable' calculates the time dynamics
%                     of the expectation value for an observable operator using the 
%                     split propagation algorithm (it is called A in the JMR paper),
%                     'trajectory' returns the system trajectory as a cell array of
%                     density matrices.
%
% The calculation of final states and observables is parallelized and tested all the 
% way to 128-core (16 nodes, 8 cores each) clusters. Parallelization of the trajectory
% calculations does not appear to yield any benefits due to large amount of inter-
% thread communication.

% Serial propagation options ('final_serial', 'observable_serial') are provided for 
% benchmarking and debugging purposes.
% 
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function answer=hs_evolution(spin_system,H,coil,rho,timestep,nsteps,output)

% Bomb out if Matlab is too old
if verLessThan('matlab', '7.12')
    error('hs_evolution: your Matlab is too old for this module - version 2011a or later required.');
end

% Disable obnoxious warning
spmd
    warning('off','distcomp:codistributed:mTimes:changeOutputCodistr');
end

% Click forward for output
spin_system=click(spin_system,'forward');

% Compute the exponential propagator
P=propagator(spin_system,H,timestep);

% Update the user
report(spin_system,'hs_evolution: propagating the system...');

% Run the evolution
switch output
    
    case 'final'

        % Parallel processing
        spmd
        
            % Slice the density matrix and the covectors
            rho=codistributed(rho,codistributor('1d',2));
            cov=speye(size(rho),codistributor('1d',2));
            
            % Propagate vectors and covectors
            for n=1:nsteps
                rho=P*rho; cov=P*cov;
            end
            
        end
        
        % Return the result
        answer=gather(rho)*gather(cov)';
        
    case 'final_dt'
        
        % Parallel processing
        spmd
            
            % Slice the density matrix column-wise
            rho=codistributed(rho,codistributor('1d',2));
            
            % Do the left side propagation
            for n=1:nsteps
                rho=P*rho;
            end
            
            % Transpose the density matrix
            rho=redistribute(rho',codistributor('1d',2));
            
            % Do the right side propagation
            for n=1:nsteps
                rho=P*rho;
            end
            
        end
        
        % Retrieve the result
        answer=gather(rho)';
        
    case 'final_serial'
        
        % Propagate the system
        for n=1:nsteps
            rho=P*rho*P';
        end
        
        % Return the result
        answer=rho;
        
    case 'observable'
        
        % Parallel processing
        spmd
        
            % Slice the density matrix and the covectors
            rho=codistributed(rho,codistributor('1d',2));
            cov=speye(size(rho),codistributor('1d',2));

            % Preallocate the local fid
            fid=zeros(nsteps+1,1);
            
            % Loop over the time steps
            for n=1:(nsteps+1)
                
                % Write the local fid
                fid(n)=sum(getLocalPart(sum(conj(cov).*(coil*rho),1)));
                
                % Step forward on rho and cov
                rho=P*rho; cov=P*cov;
                
            end
            
            % Collect the results
            answer=gplus(fid,1);
            
        end
        
        % Return the result
        answer=answer{1};
        
    case 'observable_serial'
        
        % Preallocate the answer
        answer=zeros(nsteps+1,1);
        
        % Loop over the time steps
        for n=1:(nsteps+1)
            
            % Compute the observable
            answer(n)=hadm(conj(coil),rho);
            
            % Step forward
            rho=P*rho*P';
            
        end
        
    case 'trajectory'
        
        % Preallocate the answer
        answer=cell(1,nsteps+1);
        
        % Compute the trajectory
        for n=1:(nsteps+1)
            answer{n}=rho; rho=P*rho*P';
        end
        
end

% Reenable obnoxious warning
spmd
    warning('on','distcomp:codistributed:mTimes:changeOutputCodistr');
end

end

% In 2006, Oxford's Magdalen College (where Erwin Schrodinger was a Fellow
% between 1933 and 1936) received a sum of money from a benefactor towards
% "increasing the art content of the College". A number of works were
% presented for competition, among them a beautiful stone obelisk, called
% "Monument to Knowledge" with the Schrodinger equation inscribed on it.
% The obelisk was rejected -- the inscriber had missed the bar off Planck's
% constant. As the College Governing Body put it, the equation, as written,
% "would have exploded the stone it was inscribed upon".

