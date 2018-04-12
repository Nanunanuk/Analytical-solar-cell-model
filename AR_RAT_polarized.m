function [ R,T,A ] = AR_RAT_polarized( komplex_n1,komplex_nAR,komplex_n2,wavelength,d_AR,theta,pol)
%%   Transfer Matrix Methode by Macleod
%    Input of theta in radians

if theta>pi/2
    errordlg('Input in radians','Input Error');
    exit
end 

phi1 = asin(komplex_n1 * sin(theta)/komplex_nAR);
phi2 = asin(komplex_nAR * sin(phi1)/komplex_n2);

phi1_ausgabe = phi1*180/pi;
phi2_ausgabe = phi2*180/pi;

k0 =  (2 * pi / wavelength);
delta = k0 * komplex_nAR * d_AR * cos (phi1);

y = 2.3544e-3; % [Si] "optical admittance" im vacuum
 
% "optical admittance" im Medium mit index n - % S-Polarisation (TE)
Y0_s = komplex_n1 * y * cos(theta); 
Y1_s = komplex_nAR * y * cos(phi1); 
Y2_s = komplex_n2 * y * cos(phi2); 

% "optical admittance" im Medium mit index n - % P-Polarisation (TM)
Y0_p = komplex_n1 * y * (1/cos(theta)); 
Y1_p = komplex_nAR * y * (1/cos(phi1)); 
Y2_p = komplex_n2 * y * (1/cos(phi2)); 

Matrix_s = [ cos(delta), (1i/Y1_s)*sin(delta)
               1i*Y1_s*sin(delta), cos(delta) ];
           
Matrix_p = [ cos(delta), (1i/Y1_p)*sin(delta)
               1i*Y1_p*sin(delta), cos(delta) ];
          
BC_s = Matrix_s*[1;Y2_s];    
B_s = BC_s(1);
C_s = BC_s(2);

BC_p = Matrix_p*[1;Y2_p];    
B_p = BC_p(1);
C_p = BC_p(2);

R_s = ((Y0_s * B_s - C_s)/(Y0_s * B_s + C_s)) * ((Y0_s * B_s - C_s)/(Y0_s * B_s + C_s))'; % TE - Reflexion
T_s = (4 * Y0_s * real(Y2_s) ) / ((Y0_s * B_s + C_s) * (Y0_s * B_s + C_s)');              % TE - Transmission
A_s = (4 * Y0_s * real(B_s*C_s' - Y2_s)) / ((Y0_s * B_s + C_s) * (Y0_s * B_s + C_s)');    % TE - Absorption

R_p = ((Y0_p * B_p - C_p)/(Y0_p * B_p + C_p)) * ((Y0_p * B_p - C_p)/(Y0_p * B_p + C_p))'; % TM - Reflexion
T_p = (4 * Y0_p * real(Y2_p) ) / ((Y0_p * B_p + C_p) * (Y0_p * B_p + C_p)');              % TM - Transmission
A_p = (4 * Y0_p * real(B_p*C_p' - Y2_p)) / ((Y0_p * B_p + C_p) * (Y0_p * B_p + C_p)');    % TM - Absorption

if pol=='s'
    R=R_s;
elseif pol =='p'
    R=R_p;
end

if pol=='s'
    A=A_s;
elseif pol =='p'
    A=A_p;
end

if pol=='s'
    T=T_s;
elseif pol =='p'
    T=T_p;
end


end