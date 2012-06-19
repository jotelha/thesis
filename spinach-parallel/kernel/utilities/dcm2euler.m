% Converts DCM into Euler angles (Varshalovich B convention).
%
% ilya.kuprov@oerc.ox.ac.uk

function [arg1,arg2,arg3]=dcm2euler(dcm)

%% Make it dead certain that the DCM is valid

% Grumble if the DCM is not quite orthogonal
if norm(dcm'*dcm-eye(3))>1e-6
    disp('dcm2euler: WARNING - the DCM is not orthogonal to 1e-6 tolerance.');
end

% Bomb out if the DCM is seriously non-orthogonal
if norm(dcm'*dcm-eye(3))>1e-2
    error('dcm2euler: the DCM is not orthogonal to 1e-2 tolerance, cannot proceed with conversion.');
end

% Grumble if the determinant is not quite unit
if abs(det(dcm)-1)>1e-6
    disp('dcm2euler: WARNING - DCM determinant is not unit to 1e-6 tolerance.');
end

% Bomb out if the determinant is too far from 1
if abs(det(dcm)-1)>1e-2
    error('dcm2euler: DCM determinant is not unit to 1e-2 tolerance, cannot proceed with conversion.');
end

% Warn the user if we end up close to the singularity
if (1-abs(dcm(3,3)))<1e-3
    disp('dcm2euler: WARNING - small beta angle, alpha and gamma are degenerate');
end

%% Get the Euler angles out 

% Get the beta angle out and wrap it into [0,pi]
beta=mod(acos(dcm(3,3)),pi);

% Do a brute force surface scan with respect to alpha and gamma
alphas=pi*linspace(0.05,1.95,20);
gammas=pi*linspace(0.05,1.95,20);
n_min=1; k_min=1; err_min=1;
for n=1:20
    for k=1:20
        err_current=norm(euler2dcm(alphas(n),beta,gammas(k))-dcm);
        if err_current<err_min;
            n_min=n; k_min=k; err_min=err_current;
        end
    end
end
alpha=alphas(n_min); gamma=gammas(k_min);

% Run the optimization on alpha and gamma
options=optimset('Display','off','LargeScale','off','TolX',1e-12,'TolFun',1e-12);
answer=fminunc(@(angles)norm(euler2dcm(angles(1),beta,angles(2))-dcm),[alpha gamma],options);
alpha=answer(1); gamma=answer(2);

% Wrap both angles into [0,2*pi]
alpha=mod(alpha,2*pi);
gamma=mod(gamma,2*pi);

%% Make sure the answer is valid 

% Make sure the result is good enough and bomb out if it's not
if norm(dcm-euler2dcm(alpha,beta,gamma))>1e-3
    disp(dcm); disp(euler2dcm(alpha,beta,gamma));
    error('dcm2euler: DCM to Euler conversion failed.');
end

%% Return the answer

% Adapt to the output style
if nargout==1||nargout==0
    arg1=[alpha beta gamma];
elseif nargout==3
    arg1=alpha; arg2=beta; arg3=gamma;
else
    error('dcm2euler: incorrect number of output arguments.');
end

end

% I would give the greatest sunset in the world for one sight of New York's
% skyline. Particularly when one can't see the details. Just the shapes. The
% shapes and the thought that made them. The sky over New York and the will
% of man made visible. What other religion do we need? And then people tell
% me about pilgrimages to some dank pesthole in a jungle where they go to do
% homage to a crumbling temple, to a leering stone monster with a pot belly,
% created by some leprous savage. Is it beauty and genius they want to see?
% Do they seek a sense of the sublime? Let them come to New York, stand on
% the shore of the Hudson, look and kneel. When I see the city from my
% window - no, I don't feel how small I am - but I feel that if a war came to
% threaten this, I would like to throw myself into space, over the city, and
% protect these buildings with my body.
%
% Ayn Rand, "The Fountainhead"

