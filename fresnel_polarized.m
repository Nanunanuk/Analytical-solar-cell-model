function [R, T, A]=fresnel_polarized(n1,n2,theta,pol)
%% Input:
% n1 & n2: Refractive indexes (real or complex) of the media
% Theta: Angle of incidence [rad]

cos_theta2=1/n2*sqrt(n2^2-n1^2*sin(theta)^2);
rs=(n1*cos(theta)-n2*cos_theta2)/(n1*cos(theta)+n2*cos_theta2);
rp=(n2*cos(theta)-n1*cos_theta2)/(n2*cos(theta)+n1*cos_theta2);

if pol=='s'
    R=rs*conj(rs);
elseif pol =='p'
    R=rp*conj(rp);
end

T=1-R;
A=0;
