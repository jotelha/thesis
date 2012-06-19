% Provides basic density matrix diagnostics capability -- prints the list of the
% most populated states, state vector norm, relaxation rates, cross-relaxation
% destinations etc.
%
% Ilya Kuprov               (ilya.kuprov@oerc.ox.ac.uk)
% Konstantin Pervushin      (p@trosy.com)

function state_diagnostics(spin_system,current_state,window)

% Set the defaults
if nargin==2
    window=min(20,spin_system.bas.nstates);
end

% Determine the state vector norm
current_state=full(current_state);
report(spin_system,['diagnostics: state vector norm: ' num2str(norm(current_state))]);

% Locate the most populated states and sort by coefficient amplitude
amps=abs(current_state).^2;
[amps_sorted,sorting_index]=sort(amps,1,'descend'); %#ok<ASGLU>
largest_amps=current_state(sorting_index(1:window));
largest_amp_states=full(spin_system.bas.basis(sorting_index(1:window),:));
largest_index=sorting_index(1:window);

report(spin_system,['diagnostics: ' num2str(window) ' most populated states (state, coefficient, state_number)']);
% Print at most the user-specified number of states
for n=1:window
    state_string=cell(1,spin_system.comp.nspins);
    for k=1:spin_system.comp.nspins
        switch largest_amp_states(n,k)
            case 0
                state_string{k}=' ... ';
            case 1
                state_string{k}=' T11 ';
            case 2
                state_string{k}=' T10 ';
            case 3
                state_string{k}=' T1m1';
            case 4
                state_string{k}=' T22 ';
            case 5
                state_string{k}=' T21 ';
            case 6
                state_string{k}=' T20 ';
            case 7
                state_string{k}=' T2m1';
            case 8
                state_string{k}=' T2m2';
        end
    end
    disp([cell2mat(state_string) '       ' pad(num2str(largest_amps(n),'%5.3e'),24) '  ' pad(num2str(largest_index(n)),6)]);
end

end

% One man with a dream will beat a hundred men with a grudge.

