% Converts angle-axis rotation parameters into a quaternion.
%
% gareth.charnock@oerc.ox.ac.uk

function q=aa2quat(phi,x,y,z)
	
    % Normalize the axis ort
    norm = sqrt(x*x+y*y+z*z);
	x=x/norm; y=y/norm;	z=z/norm;

    % Compute the quaternion
	q.w=cos(a/2);
	q.x=x*sin(phi/2);
	q.y=y*sin(phi/2);
	q.z=z*sin(phi/2);
    
end

% You'll get everything society can give a man. You'll keep all the money.
% You'll take any fame or honor anyone might want to grant. You'll accept
% such gratitude as the tenants might feel. And I - I'll take what nobody
% can give a man, except himself. I will have built Cortlandt.
%
% Ayn Rand, "The Fountainhead"

