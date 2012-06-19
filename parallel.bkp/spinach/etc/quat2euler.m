% Converts Quaternion into Euler angles (Varshalovich B convention).
%
% Quaternion is defined as a structure with the following fields: q.w,q.x,q.y,q.z
%
% gareth.charnock@oerc.ox.ac.uk

function [alpha,beta,gamma]=quat2euler(q)

    % Normalise the quaternion
    norm=sqrt(q.x^2+q.y^2+q.z^2+q.w^2);
    q.x=q.x/norm;
    q.y=q.y/norm;
    q.z=q.z/norm;
    q.w=q.w/norm;
    
    %Peform the conversion
    z_axis=[0 0 1];
    x_axis=[1 0 0];
    
    z_axis=quaternion_transform(q,z_axis);
    x_axis=quaternion_transform(q,x_axis);
    
    gamma=atan2(z_axis(2),z_axis(1));
    beta=atan2(sqrt(z_axis(1)^2 + z_axis(2)^2),z_axis(3));

    betaTwist=aa2quat(-beta,0,1,0);
    gammaTwist=aa2quat(-gamma,0,0,1);

    x_axis=quaternion_transform(mult_quaternion(betaTwist,gammaTwist),x_axis);
    alpha = atan2(x_axis(2),x_axis(1));

    alpha = mod(alpha,2*pi);
    beta  = mod(beta,   pi);
    gamma = mod(gamma,2*pi);
    
end

function q_star=quaternion_conjugulate(q)
    q_star.w=q.w;
    q_star.x=-q.x;
    q_star.y=-q.y;
    q_star.z=-q.z;
end

function v_prime = quaternion_transform(q,v)
    q_vec.w = 0;
    q_vec.x = v(1);
    q_vec.y = v(2);
    q_vec.z = v(3);
    
    q_vec=mult_quaternion(q,mult_quaternion(q_vec,quaternion_conjugulate(q)));
    
    v_prime(1) = q_vec.x;
    v_prime(2) = q_vec.y;
    v_prime(3) = q_vec.z;
end

function [q]=mult_quaternion(q1,q2)
	q.w=q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
	q.x=q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
	q.y=q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x;
	q.z=q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w;	
end

% The paper describing the time-resolved photo-CIDNP installation at Oxford
% (http://link.aip.org/link/doi/10.1063/1.2010287) was written by IK whilst
% under an indefinite ban from using said installation. The ban was put in
% place after Keith McLauchlan wandered idly into the NMR lab one day and
% discovered a running Nd:YAG laser (18 megawatt pulse amplitude, UV) tucked
% onto a windowsill and three quartz prisms bolted at odd angles to walls
% and ceiling, directing the beam into the NMR magnet from above. KMcL bolted 
% from the lab in the direction of the Safety Officer and, in ten minutes,
% visibly shaken Gus Hancock pronounced the lab closed for operation until
% the beam of "that Russian contraption" was properly enclosed. Thankfully,
% enough data was collected at that point for three more papers - the system
% was subsequently dismantled as surplus to requirement.

