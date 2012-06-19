% GRAPE objective function and gradient. Propagates the system through a 
% user-supplied shaped pulse from a given initial state and projects the
% result onto the given final state. The real part of the projection is 
% returned, along with its gradient with respect to amplitudes of all ope-
% rators in every step of the shaped pulse. The input parameters are:
%
%      drift          -       the "drift" Liouvillian, i.e. the couplings,
%                             relaxation and other such things that continue
%                             operating while the pulse is being executed.
%
%      controls       -       a cell array of commutation superoperators whose
%                             coefficients are getting modulated by the pulse,
%                             e.g. {Lx,Ly}
%
%      waveform       -       matrix of doubles with n_controls rows and
%                             n_steps columns, where n_controls is the number
%                             of control operators and n_steps is the number
%                             of steps in the waveform.
%
%      time_step      -       waveform time step, in seconds.
%
%      nsteps         -       number of steps in the waveform.
%
%      starting_state_array - a cell array of initial state vectors (if
%                             multiple sources and targets are supplied, the
%                             resulting objective function and the gradient
%                             are sums over these source-target pairs).
%
%      target_state_array   - a cell array of target state vectors (if multi-
%                             ple sources and targets are supplied, the resul-
%                             ting objective function and the gradient are 
%                             sums over these source-target pairs).
%
%      scaling_factor       - internal scaling factor. Set this to match
%                             the average power level, this improves the 
%                             numerical acuracy and optimizer performance.
%
% Three derivative calculation methods are available: Hausdorff series, second
% order central funite difference and fourth order central finite difference.
% The default (second order CFD) is sufficiently accurate for most purposes.
%                             
% ilya.kuprov@oerc.ox.ac.uk
%
%#ok<*PFBNS>

function [total_objective,total_grad]=grape(spin_system,drift,controls,waveform,time_step,nsteps,starting_state_array,target_state_array,scaling_factor)

% Set the default scaling factor
if ~exist('scaling_factor','var')
    scaling_factor=1;
end

% Allow numeric inputs for single-source, single-target calculations
if isnumeric(starting_state_array), starting_state_array={starting_state_array}; end
if isnumeric(target_state_array), target_state_array={target_state_array}; end

% Preallocate the results
total_objective=0; total_grad=zeros(size(waveform));

% Scale the waveform
waveform=scaling_factor*waveform;

% Loop over the sources and targets
for contribution=1:numel(starting_state_array)
    
    % Grab the current source and target
    starting_state=starting_state_array{contribution};
    target_state=target_state_array{contribution};
    
    % Preallocate the trajectory arrays
    forward_traj=zeros(size(starting_state,1),nsteps+1);
    backward_traj=zeros(size(starting_state,1),nsteps+1);
    
    % Run the forward propagation
    forward_traj(:,1)=starting_state;
    for n=1:nsteps
        L=drift;
        for k=1:length(controls)
            L=L+waveform(k,n)*controls{k};
        end
        forward_traj(:,n+1)=step(spin_system,L,forward_traj(:,n),time_step);
    end
    
    % Run the backward propagation
    backward_traj(:,nsteps+1)=target_state;
    for n=nsteps:-1:1
        L=drift';
        for k=1:length(controls)
            L=L+waveform(k,n)*controls{k};
        end
        backward_traj(:,n)=step(spin_system,L,backward_traj(:,n+1),-time_step);
    end
    
    % Compute the gradient
    if nargout==2
        
        % Preallocate the array
        grad=zeros(size(waveform));
        
        % Loop over time steps
        parfor n=1:nsteps
            
            % Initialize the Liouvillian
            L=drift;
            
            % Add the current controls
            for k=1:length(controls)
                L=L+waveform(k,n)*controls{k};
            end
            A=-1i*L*time_step;
            
            % Preallocate the derivatives
            control_derivatives=zeros(length(controls),1);
            
            % Generate the control derivatives
            switch spin_system.tols.dP_method
                
                case 'hausdorff'
                    
                    % Detect the abuse
                    if norm(A,'inf')>10
                        error('grape: waveform time step greatly exceeds the scale of system dynamics.');
                    end
                    
                    % Loop over the controls
                    for k=1:length(controls)
                        
                        % Converge the Hausdorff series
                        H=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,spin_system.bas.nstates);
                        next_term=-1i*controls{k}*time_step; m=1;
                        while norm(next_term,'inf') > spin_system.tols.derivative_prop_norm
                            
                            % Add to the total
                            H=H+next_term; m=m+1;
                            
                            % Obey the accuracy target
                            if m > spin_system.tols.dP_order, break; end
                            
                            % Compute the next term
                            next_term=(next_term*A-A*next_term)/m;
                            
                            % Run the clean-up procedure and check the density
                            next_term=clean_up(spin_system,next_term,spin_system.tols.derivative_prop_chop);
                            
                        end
                        
                        % Compute the control derivatives
                        control_derivatives(k)=backward_traj(:,n)'*H*forward_traj(:,n);
                        
                    end
                    
                case 'fd_O(h^2)'
                    
                    % Get the finite difference step
                    delta=max(abs(waveform(:)))*sqrt(eps);
                    
                    % Loop over the controls
                    for k=1:length(controls)
                        
                        % Compute the control derivatives with central finite differences
                        control_derivatives(k)=backward_traj(:,n+1)'*(step(spin_system,L+delta*controls{k},forward_traj(:,n),time_step)-...
                                                                      step(spin_system,L-delta*controls{k},forward_traj(:,n),time_step))/(2*delta);
                        
                    end
                    
                case 'fd_O(h^4)'
                    
                    % Get the finite difference step
                    delta=max(abs(waveform(:)))*sqrt(eps);
                    
                    % Loop over the controls
                    for k=1:length(controls)
                        
                        % Compute the control derivatives with central finite differences
                        control_derivatives(k)=backward_traj(:,n+1)'*(-1*step(spin_system,L+2*delta*controls{k},forward_traj(:,n),time_step)...
                                                                      +8*step(spin_system,L+delta*controls{k},forward_traj(:,n),time_step)...
                                                                      -8*step(spin_system,L-delta*controls{k},forward_traj(:,n),time_step)...
                                                                      +1*step(spin_system,L-2*delta*controls{k},forward_traj(:,n),time_step))/(12*delta);
                        
                    end
                    
                otherwise
                    
                    error('grape: unknown propagator differentiation method.');
                    
            end
            
            % Assign the gradients for the current step
            grad(:,n)=control_derivatives;
            
        end
        
        % Add the gradient to the total
        total_grad=total_grad+real(grad);
        
    end
    
    % Compute the cost function and add to the total
    total_objective=total_objective+real(target_state'*forward_traj(:,end));
    
end

% Normalize the total objective and the total gradient
total_objective=total_objective/numel(starting_state_array);
total_grad=scaling_factor*total_grad/numel(starting_state_array);

end

% In any culture, subculture, or family in which belief is valued above
% thought, self-surrender is valued above self-expression, and conformity
% is valued above integrity, those who preserve their self-esteem are likely 
% to be heroic exceptions.
%
% Nathaniel Branden

