clear all
clc

store_G=0;

%% User has to select a folder where the generation profiles for PC1D are stored 
if store_G==1
    [save_file_loc] = uigetdir('C:\Users\Nils\Documents\MATLAB\Grundlagen\Reflection\Reflection_including_R_esc\Results','Select path to store output data');
end

%% Load material parameters

[a b c]=xlsread(['C:\Users\Nils\Documents\Nils Reiners\Messdaten\Material_And_Spectral_Data\PV-Lighthouse' '\All_data_PV-Lighthouse.xlsx']);

Mat(1,:)={'Air' 'Air' '[Pal85a]'};
Mat(2,:)={'Si' 'Crystalline 300 K' '[Gre08]'};
Mat(3,:)={'Al' 'Pure' '[Pal85a]'};
Mat(4,:)={'Ag','Pure','[Pal85a]'};
Mat(5,:)={'SiNx','PECVD','[Bak11]'};
Mat(6,:)={'SiO2','Thermal','[Pal85e]'};
Mat(7,:)={'AlSi', 'Eutectic', '[Vog15]'};
Mat(8,:)={'Glass','Low-Fe Starphire','[McI09b]'};

for n = 1:size(Mat,1)
    row1=find(strcmp(Mat{n,1},b(1,:)));
    row2=find(strcmp(Mat{n,2},b(2,:)));
    row3=find(strcmp(Mat{n,3},b(3,:)));
    
    x1=row1(ismember(row1,row2));
    data_col=x1(ismember(x1,row3));
    
    wl_n_k=[a(:,1) a(:,data_col) a(:,data_col+1)];
    assignin('base', [Mat{n}], wl_n_k);
end


%% Read in solar spectrum
[a b c]=xlsread(['C:\Users\Nils\Documents\Nils Reiners\Messdaten\Material_And_Spectral_Data\PV-Lighthouse' '\PVL_Spectrum_AM15g.xlsx'],'Spectrum');
AM15=[a(:,1) a(:,2)];
%% Basic constants
epsilon=11.8*8.854187817e-12;   %[As/Vm]
q=1.602E-19;                    %[Js]
k=1.38064852E-23;               %[J/K]
T=300;                          %[K]
n_i=1e10;
h=6.62606957E-34;
c_0=299792458;
%% Global parameter


batch=1;
batch_init=1;

for theta_0=[8 70];                    %[�]
    
    %% Creating folder for PC1D generation files
    if store_G==1
        formatOut = 'yyyy_mm_dd_hh_MM';
        folderName=[datestr(now,formatOut) '_gen_batch-' num2str(batch) '_angle-' num2str(theta_0)];
        mkdir(save_file_loc,folderName)
        Gen_folder=[save_file_loc '\' folderName];
    end
    
    d_AR=87;                        %[nm]
    for W=[100 200]./1000000   ;               %[m]
        lambda=(300:10:1200)';          %[nm]
        
        % Lambertian share and reflection intensity
        for d_rel_b=[0];     
            for d_rel_f=[0];
                for r_db=[0.8];
                    for r_df=[0.98];
                        
                        % Grid
                        pf=2150;
                        
                        wf= 0;%180;                                                % If effective fingerwidth is decreasing with AOI: 100*1/(cos(pi/3*theta_0/90));

                        
                        pbb=60.07; %[mm]
                        wbb=0;
                        
                        Mbb=wbb/(pbb+wbb);
                        Mfi=pbb*wf/((pbb+wbb)*(pf+wf));
                        Msi=pbb*pf/((pbb+wbb)*(pf+wf));
                        
                        N_D=1e20;
                        N_A=1e16;
                        
                        V_0= T*k/q*log(N_A*N_D/n_i^2);
                        
                        for W_e=[0.04] /1000000;                               %[m]
                            W_scr=sqrt(2*epsilon/q*V_0*(1/N_A+1/N_D))/1000; %[m]
                            W_b=W-W_e-W_scr;                                %[m]                     
                            
                            L_n_count=1;
                            for L_n=[100 200]   /1000000;
                            L_n_contour(L_n_count)=L_n*1000000;
                            L_n_count=L_n_count+1;
                            
                                for L_p=[0.4]     /1000000;
                                    for S_e=[2E5] /100;
                                        
                                        S_b_count=1;
                                        for S_b=[0 10000] /100;
                                            S_b_contour(S_b_count)=S_b*100;
                                            S_b_count=S_b_count+1;
                                
                                            D_p=5   /10000;
                                            D_n=29  /10000;
                                            
                                            alpha=4*pi*interp1(Si(:,1),Si(:,3),lambda)./(lambda(:,1)/1e9);
                                            phi=interp1(AM15(:,1),AM15(:,2),lambda)/(h*c_0*1E9).*lambda;
                                            
                                            %% Calculation for s and p polarized light
                                            
                                            for pol=['s' 'p'];
                                                
                                                %% Starting the simulation for all Wavelengths
                                                for n=1:length(lambda)
                                                    %% Free carrier absorption
                                                    
                                                    alpha_FCA_e(n,:)=0.00000168*N_D*(lambda(n)/10000000)^2.88*100;
                                                    alpha_FCA_scr(n,:)=0.00000168*(1/2*(N_D+n_i))*(lambda(n)/10000000)^2.88*100;
                                                    alpha_FCA_b(n,:)=0.0000000018*N_A*(lambda(n)/10000000)^2.18*100;
                                                    
                                                    p_si(n,:)=alpha(n,:)*W...
                                                        /(alpha(n,:)*W...
                                                        +alpha_FCA_e(n,:)*W_e...
                                                        +alpha_FCA_scr(n,:)*W_scr...
                                                        +alpha_FCA_b(n,:)*W_b);
                                                    
                                                    p_FCA(n,:)=(alpha_FCA_e(n,:)*W_e...
                                                        +alpha_FCA_scr(n,:)*W_scr...
                                                        +alpha_FCA_b(n,:)*W_b)...
                                                        /(alpha(n,:)*W...
                                                        +alpha_FCA_e(n,:)*W_e...
                                                        +alpha_FCA_scr(n,:)*W_scr...
                                                        +alpha_FCA_b(n,:)*W_b);
                                                    
                                                    %% Internal angles
                                                    
                                                    theta_1(n,:)= asin(sin(theta_0*pi/180)*Air(n,2)/Si(n,2));        %[rad]  %For 45� texture emulation (10/45*theta_0+35)*pi/180;

                                                    theta_2(n,:)=theta_1(n,:);                                      %[rad]
                                                    theta_n(n,:)=60*pi/180;                                         %[rad]
                                                    
                                                    %% Testing for angular dependent d_rel for front and backside
                                                    
                                                    %d_rel_b=(exp(theta_1(n,:)*log(d_rel_b_0+1)/(pi/2))-1);%+d_rel_b_0;
                                                    
                                                    
                                                    %d_rel_f=-(exp(theta_1(n,:)*log(d_rel_f_0+1)/(pi/2))-1)+d_rel_f_0;
                                                    %% RAT_fe
                                                    
                                                    n_1=interp1(Air(:,1),Air(:,2),lambda(n));
                                                    komplex_n_AR = interp1(SiNx(:,1),SiNx(:,2 )- SiNx(:,3) * 1i,lambda(n));
                                                    komplex_n_2 = interp1(Si(:,1),Si(:,2) - Si(:,3) * 1i,lambda(n));
                                                    
                                                    [r_fe(n,:), t_fe(n,:), a_fe(n,:)] = AR_RAT_polarized(n_1,komplex_n_AR,komplex_n_2,lambda(n),d_AR,theta_0*pi/180,pol);
                                                    
                                                    %% RAT_b1
                                                    n_1=interp1(Si(:,1),Si(:,2),lambda(n));
                                                    n_2 = interp1(SiO2(:,1),SiO2(:,2),lambda(n))-interp1(SiO2(:,1),SiO2(:,3),lambda(n))*i;
                                                    n_3 = interp1(Al(:,1),Al(:,2),lambda(n))-interp1(Al(:,1),Al(:,3),lambda(n))*i;
                                                    d_AlSi=5;
                                                    
                                                    [r_b1(n,:), t_b1(n,:), a_b1(n,:)] = AR_RAT_polarized(n_1,n_2,n_3,lambda(n),d_AlSi,theta_1(n),pol);
                                                    
                                                    
                                                    %% RAT_f1
                                                    n_1= interp1(Si(:,1),Si(:,2) - Si(:,3) * 1i,lambda(n));
                                                    komplex_n_AR = interp1(SiNx(:,1),SiNx(:,2 )- SiNx(:,3) * 1i,lambda(n));
                                                    komplex_n_2 = interp1(Air(:,1),Air(:,2),lambda(n));
                                                    
                                                    [r_f1(n,:), t_f1(n,:), a_f1(n,:)] = AR_RAT_polarized(n_1,komplex_n_AR,komplex_n_2,lambda(n),d_AR,theta_1(n),pol);
                                                    
                                                    %% RAT_fi
                                                    n_1 = 1;
                                                    n_2 = interp1(Ag(:,1),Ag(:,2),lambda(n))-interp1(Ag(:,1),Ag(:,3),lambda(n))*i;
                                                    
                                                    [r_fi(n,:), t_fi(n,:), a_fi(n,:)] = fresnel_polarized(n_1,n_2,theta_0*pi/180,pol);
                                                    
                                                    
                                                    %% Calculation of major RAT values
                                                    R_fe(n,:)=Msi*r_fe(n,:)+Mfi*r_fi(n,:)+Mbb;
                                                    
                                                    T_fe(n,:)=Msi*t_fe(n,:);
                                                    
                                                    A_grid(n,:)=Mfi*(1-r_fi(n,:));
                                                    
                                                    A_AR(n,:)=Msi*a_fe(n,:);
                                                    
                                                    if theta_0==8
                                                        R_b1(n,:)=r_b1(n,:);
                                                    else
                                                        R_b1(n,:)=r_b1(n,:);
                                                    end
                                                    
                                                    if theta_0==8
                                                        R_f1(n,:)=0;%Mfi+Msi*r_f1(n,:);
                                                    else
                                                        R_f1(n,:)=0;%Mfi+Msi*r_f1(n,:);
                                                    end
                                                    
                                                    T_f1(n,:)=1-R_f1(n,:);
                                                    
                                                    
                                                    %% Calculation of direct transmission through the wafer
                                                    T_d1(n,:)=exp(-(alpha(n,:)+alpha_FCA_e(n,:))*W_e/cos(theta_1(n,:)))...
                                                        *exp(-(alpha(n,:)+alpha_FCA_scr(n,:))*W_scr/cos(theta_1(n,:)))...
                                                        *exp(-(alpha(n,:)+alpha_FCA_b(n,:))*W_b/cos(theta_1(n,:)));
                                                    
                                                    
                                                    %% Diffuse transmission through the wafer with angle of 60�
                                                    T_diff(n,:)=exp(-(alpha(n,:)+alpha_FCA_e(n,:))*W_e/cos(60*pi/180))...
                                                        *exp(-(alpha(n,:)+alpha_FCA_scr(n,:))*W_scr/cos(60*pi/180))...
                                                        *exp(-(alpha(n,:)+alpha_FCA_b(n,:))*W_b/cos(60*pi/180));
                                                    
                                                    %% Reflection ("realistic diffuse reflection model")
                                                    R_esc_dir(n,:)=T_fe(n,:)*(1-d_rel_b)*(1-d_rel_f)*T_d1(n,:)^2*R_b1(n,:)*T_f1(n,:)...
                                                        /(1-(1-d_rel_b)*(1-d_rel_f)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:));
                                                    
                                                    R_esc_diff_f(n,:)=T_fe(n,:)*d_rel_f*T_diff(n,:)^2*r_db*(1-r_df)...
                                                        /((T_diff(n,:)^2*r_db*r_df-1)*((1-d_rel_b)*(1-d_rel_f)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:)-1));
                                                    
                                                    R_esc_diff_b(n,:)=T_fe(n,:)*d_rel_b*(1-d_rel_f)*T_diff(n,:)*T_d1(n,:)*R_b1(n,:)*(1-r_df)...
                                                        /((T_diff(n,:)^2*r_db*r_df-1)*((1-d_rel_b)*(1-d_rel_f)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:)-1));
                                                    
                                                    R_esc(n,:)=R_esc_dir(n,:)+R_esc_diff_f(n,:)+R_esc_diff_b(n,:);
                                                    
                                                    R_tot(n,:)=R_fe(n,:)+R_esc(n,:);
                                                    
                                                    %% Absorption ("realistic diffuse reflection model")
                                                    
                                                    A_dir(n,:)=T_fe(n,:)*(1-d_rel_f)*(1+(1-d_rel_b)*T_d1(n,:)*R_b1(n,:))...
                                                        *(1-T_d1(n,:))...
                                                        /(1-(1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:));
                                                    
                                                    A_diff(n,:)=T_fe(n,:)...
                                                        *(d_rel_f...
                                                        +d_rel_b*(1-d_rel_f)*T_diff(n,:)*T_d1(n,:)*r_df*R_b1(n,:)...
                                                        +d_rel_f*T_diff(n,:)*r_db...
                                                        +d_rel_b*(1-d_rel_f)*T_d1(n,:)*R_b1(n,:))...
                                                        *(1-T_diff(n,:))...
                                                        /((T_diff(n,:)^2*r_db*r_df-1)*((1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:)-1));
                                                    
                                                    A_tot(n,:)=A_dir(n,:)+A_diff(n,:);
                                                    
                                                    A_si(n,:)=p_si(n,:)...
                                                        *A_tot(n,:);
                                                    
                                                    A_FCA(n,:)=p_FCA(n,:)...
                                                        *A_tot(n,:);
                                                    
                                                    %% Absorption in the backside = Transmission ("realistic diffuse reflection model")
                                                    
                                                    A_back_dir(n,:)=T_fe(n,:)*(1-d_rel_f)*T_d1(n,:)*(1-R_b1(n,:))...
                                                        /(1-(1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:));
                                                    
                                                    A_back_diff(n,:)=T_fe(n,:)*(d_rel_f*T_diff(n,:)+d_rel_b*(1-d_rel_f)*T_diff(n,:)^2*T_d1(n,:)*r_df*R_b1(n,:))*(1-r_db)...
                                                        /((T_diff(n,:)^2*r_db*r_df-1)*((1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:)-1));
                                                    
                                                    A_back_tot(n,:)=A_back_dir(n,:)+A_back_diff(n,:);
                                                    
                                                    %% Cumulative photon flux
                                                    
                                                    phi_hat_a(n,:)=phi(n,:)*T_fe(n,:)*(1-d_rel_f)...
                                                        /(1-(1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:));
                                                    phi_hat_b(n,:)=phi(n,:)*T_fe(n,:)*(1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)*R_b1(n,:)...
                                                        /(1-(1-d_rel_f)*(1-d_rel_b)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:));
                                                    phi_hat_c(n,:)=phi(n,:)*T_fe(n,:)...
                                                        *(d_rel_f+d_rel_b*(1-d_rel_f)*T_diff(n,:)*T_d1(n,:)*r_df*R_b1(n,:))...
                                                        /((T_diff(n,:)^2*r_db*r_df-1)...
                                                        *((1-d_rel_b)*(1-d_rel_f)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:)-1));
                                                    phi_hat_d(n,:)=phi(n,:)*T_fe(n,:)...
                                                        *(d_rel_f*T_diff(n,:)*r_db+d_rel_b*(1-d_rel_f)*T_d1(n,:)*R_b1(n,:))...
                                                        /((T_diff(n,:)^2*r_db*r_df-1)...
                                                        *((1-d_rel_b)*(1-d_rel_f)*T_d1(n,:)^2*R_b1(n,:)*R_f1(n,:)-1));
                                                    
                                                    %% Monochromatic generation profile for the whole cell
                                                    
                                                    if store_G==1
                                                        x1=depth_array(0.001,W*1000000,10000,300)/1000000; %Creating an array of depth values with shorter distances close to the surfaces, units=[m]
                                                        
                                                        % The following model is using the phi hat equations from above
                                                        G_cell(n,:)=alpha(n,:)/cos(theta_1(n,:))...
                                                            *p_si(n,:)...
                                                            *(phi_hat_a(n,:)*exp(-alpha(n,:)/cos(theta_1(n,:))*x1)...
                                                            +phi_hat_b(n,:)*exp(-alpha(n,:)/cos(theta_1(n,:))*(W-x1)))...
                                                            +2*alpha(n,:)*(phi_hat_c(n,:)*exp(-2*alpha(n,:)*x1)...
                                                            +phi_hat_d(n,:)*exp(-2*alpha(n,:)*(W-x1)));
                                                        
                                                        %% Creating cummulative generation profile
                                                    
                                                        for ii=1:length(G_cell(n,:))
                                                            if ii==1
                                                                G_cum_long(1,:)=zeros(1,length(G_cell(n,:)));
                                                            else
                                                                G_cum_long(ii)=(x1(ii)-x1(ii-1))*G_cell(n,ii)+G_cum_long(ii-1);
                                                            end
                                                        end
                                                        x2=depth_array(0.001,W*1000000,90,lambda(n))/1000000;
                                                        G_cum=interp1(x1,G_cum_long,x2);
                                                        assignin('base', ['G_cum_' num2str(lambda(n))], [x2'*1000000 G_cum']);
                                                        
                                                        A_test(n)=max(G_cum);
                                                        
                                                        %% Creating input file for PC1D
                                                        
                                                        cd(Gen_folder);
                                                        fid = fopen(['Gen_' num2str(lambda(n)) '.gen'], 'wt');
                                                        for ii=1:length(G_cum)
                                                            fprintf(fid, '%f\t' , [x2(ii)*1000000 G_cum(ii)/100]);
                                                            fprintf(fid, '\n');
                                                        end
                                                        fclose(fid);
                                                    end
                                                    %% Generation spectrum inside the emitter
                                                    f_ae(n,:)=phi_hat_a(n,:)...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_e(n,:))/cos(theta_1(n,:))*W_e));
                                                    f_be(n,:)=phi_hat_b(n,:)...
                                                        *exp(-((alpha(n,:)+alpha_FCA_b(n,:))*W_b+(alpha(n,:)+alpha_FCA_scr(n,:))*W_scr)/cos(theta_2(n,:)))...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_e(n,:))/cos(theta_2(n,:))*W_e));
                                                    f_ce(n,:)=phi_hat_c(n,:)...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_e(n,:))/cos(theta_n(n,:))*W_e));
                                                    f_de(n,:)=phi_hat_d(n,:)...
                                                        *exp(-((alpha(n,:)+alpha_FCA_b(n,:))*W_b+(alpha(n,:)+alpha_FCA_scr(n,:))*W_scr)/cos(theta_n(n,:)))...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_e(n,:))/cos(theta_n(n,:))*W_e));
                                                    
                                                    %% Generation spectrum inside the scr
                                                    f_ascr(n,:)=phi_hat_a(n,:)...
                                                        *exp(-(alpha(n,:)+alpha_FCA_e(n,:))/cos(theta_1(n,:))*W_e)...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_scr(n,:))/cos(theta_1(n,:))*W_scr));
                                                    f_bscr(n,:)=phi_hat_b(n,:)...
                                                        *exp(-(alpha(n,:)+alpha_FCA_b(n,:))*W_b/cos(theta_2(n,:)))...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_scr(n,:))/cos(theta_2(n,:))*W_scr));
                                                    f_cscr(n,:)=phi_hat_c(n,:)...
                                                        *exp(-(alpha(n,:)+alpha_FCA_e(n,:))/cos(theta_n(n,:))*W_e)...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_scr(n,:))/cos(theta_n(n,:))*W_scr));
                                                    f_dscr(n,:)=phi_hat_d(n,:)...
                                                        *exp(-(alpha(n,:)+alpha_FCA_b(n,:))*W_b/cos(theta_n(n,:)))...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_scr(n,:))/cos(theta_n(n,:))*W_scr));
                                                    
                                                    %% Generation spectrum inside the base
                                                    f_ab(n,:)=phi_hat_a(n,:)...
                                                        *exp(-((alpha(n,:)+alpha_FCA_e(n,:))*W_e+(alpha(n,:)+alpha_FCA_b(n,:))*W_scr)/cos(theta_1(n,:)))...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_b(n,:))/cos(theta_1(n,:))*W_b));
                                                    f_bb(n,:)=phi_hat_b(n,:)...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_b(n,:))/cos(theta_2(n,:))*W_b));
                                                    f_cb(n,:)=phi_hat_c(n,:)...
                                                        *exp(-((alpha(n,:)+alpha_FCA_b(n,:))*W_e+(alpha(n,:)+alpha_FCA_b(n,:))*W_scr)/cos(theta_n(n,:)))...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_b(n,:))/cos(theta_n(n,:))*W_b));
                                                    f_db(n,:)=phi_hat_d(n,:)...
                                                        *p_si(n,:)...
                                                        *(1-exp(-(alpha(n,:)+alpha_FCA_b(n,:))/cos(theta_n(n,:))*W_b));
                                                    
                                                    %% Spectrum of the collection probability inside the emitter
                                                    eta_forw_e_a(n,:)=alpha(n,:)/cos(theta_1(n,:))*L_p/...
                                                        (((alpha(n,:)/cos(theta_1(n,:)))^2*L_p^2-1)*(1-exp(-alpha(n,:)/cos(theta_1(n,:))*W_e)))...
                                                        *(alpha(n,:)/cos(theta_1(n,:))*L_p...
                                                        -(sinh(W_e/L_p)+S_e*L_p/D_p*cosh(W_e/L_p)-(S_e*L_p/D_p-alpha(n,:)/cos(theta_1(n,:))*L_p)...
                                                        *exp(-alpha(n,:)/cos(theta_1(n,:))*W_e))...
                                                        /(cosh(W_e/L_p)+S_e*L_p/D_p*sinh(W_e/L_p)));
                                                    eta_back_e_b(n,:)=alpha(n,:)/cos(theta_2(n,:))*L_p/...
                                                        (((alpha(n,:)/cos(theta_2(n,:)))^2*L_p^2-1)*(1-exp(-alpha(n,:)/cos(theta_2(n,:))*W_e)))...
                                                        *(-alpha(n,:)/cos(theta_2(n,:))*L_p*exp(-alpha(n,:)/cos(theta_2(n,:))*W_e)-...
                                                        ((sinh(W_e/L_p)+S_e*L_p/D_p*cosh(W_e/L_p))*exp(-alpha(n,:)/cos(theta_2(n,:))*W_e)...
                                                        -(S_e*L_p/D_p+alpha(n,:)/cos(theta_2(n,:))*L_p))...
                                                        /(cosh(W_e/L_p)+S_e*L_p/D_p*sinh(W_e/L_p)));
                                                    eta_forw_e_c(n,:)=alpha(n,:)/cos(theta_n(n,:))*L_p/...
                                                        (((alpha(n,:)/cos(theta_n(n,:)))^2*L_p^2-1)*(1-exp(-alpha(n,:)/cos(theta_n(n,:))*W_e)))...
                                                        *(alpha(n,:)/cos(theta_n(n,:))*L_p-...
                                                        (sinh(W_e/L_p)+S_e*L_p/D_p*cosh(W_e/L_p)-(S_e*L_p/D_p-alpha(n,:)/cos(theta_n(n,:))*L_p)...
                                                        *exp(-alpha(n,:)/cos(theta_n(n,:))*W_e))...
                                                        /(cosh(W_e/L_p)+S_e*L_p/D_p*sinh(W_e/L_p)));
                                                    eta_back_e_d(n,:)=alpha(n,:)/cos(theta_n(n,:))*L_p/...
                                                        (((alpha(n,:)/cos(theta_n(n,:)))^2*L_p^2-1)*(1-exp(-alpha(n,:)/cos(theta_n(n,:))*W_e)))...
                                                        *(-alpha(n,:)/cos(theta_n(n,:))*L_p*exp(-alpha(n,:)/cos(theta_n(n,:))*W_e)-...
                                                        ((sinh(W_e/L_p)+S_e*L_p/D_p*cosh(W_e/L_p))*exp(-alpha(n,:)/cos(theta_n(n,:))*W_e)...
                                                        -(S_e*L_p/D_p+alpha(n,:)/cos(theta_n(n,:))*L_p))...
                                                        /(cosh(W_e/L_p)+S_e*L_p/D_p*sinh(W_e/L_p)));
                                                    
                                                    %% Spectrum of the collection probability inside the base
                                                    eta_forw_b_a(n,:)=alpha(n,:)/cos(theta_1(n,:))*L_n/...
                                                        (((alpha(n,:)/cos(theta_1(n,:)))^2*L_n^2-1)*(1-exp(-alpha(n,:)/cos(theta_1(n,:))*W_b)))...
                                                        *(alpha(n,:)/cos(theta_1(n,:))*L_n-...
                                                        (sinh(W_b/L_n)+S_b*L_n/D_n*cosh(W_b/L_n)-(S_b*L_n/D_n-alpha(n,:)/cos(theta_1(n,:))*L_n)...
                                                        *exp(-alpha(n,:)/cos(theta_1(n,:))*W_b))...
                                                        /(cosh(W_b/L_n)+S_b*L_n/D_n*sinh(W_b/L_n)));
                                                    eta_back_b_b(n,:)=alpha(n,:)/cos(theta_2(n,:))*L_n/...
                                                        (((alpha(n,:)/cos(theta_2(n,:)))^2*L_n^2-1)*(1-exp(-alpha(n,:)/cos(theta_2(n,:))*W_b)))...
                                                        *(-alpha(n,:)/cos(theta_2(n,:))*L_n*exp(-alpha(n,:)/cos(theta_2(n,:))*W_b)-...
                                                        ((sinh(W_b/L_n)+S_b*L_n/D_n*cosh(W_b/L_n))*exp(-alpha(n,:)/cos(theta_2(n,:))*W_b)...
                                                        -(S_b*L_n/D_n+alpha(n,:)/cos(theta_2(n,:))*L_n))...
                                                        /(cosh(W_b/L_n)+S_b*L_n/D_n*sinh(W_b/L_n)));
                                                    eta_forw_b_c(n,:)=alpha(n,:)/cos(theta_n(n,:))*L_n/...
                                                        (((alpha(n,:)/cos(theta_n(n,:)))^2*L_n^2-1)*(1-exp(-alpha(n,:)/cos(theta_n(n,:))*W_b)))...
                                                        *(alpha(n,:)/cos(theta_n(n,:))*L_n-...
                                                        (sinh(W_b/L_n)+S_b*L_n/D_n*cosh(W_b/L_n)-(S_b*L_n/D_n-alpha(n,:)/cos(theta_n(n,:))*L_n)...
                                                        *exp(-alpha(n,:)/cos(theta_n(n,:))*W_b))...
                                                        /(cosh(W_b/L_n)+S_b*L_n/D_n*sinh(W_b/L_n)));
                                                    eta_back_b_d(n,:)=alpha(n,:)/cos(theta_n(n,:))*L_n/...
                                                        (((alpha(n,:)/cos(theta_n(n,:)))^2*L_n^2-1)*(1-exp(-alpha(n,:)/cos(theta_n(n,:))*W_b)))...
                                                        *(-alpha(n,:)/cos(theta_n(n,:))*L_n*exp(-alpha(n,:)/cos(theta_n(n,:))*W_b)-...
                                                        ((sinh(W_b/L_n)+S_b*L_n/D_n*cosh(W_b/L_n))*exp(-alpha(n,:)/cos(theta_n(n,:))*W_b)...
                                                        -(S_b*L_n/D_n+alpha(n,:)/cos(theta_n(n,:))*L_n))...
                                                        /(cosh(W_b/L_n)+S_b*L_n/D_n*sinh(W_b/L_n)));
                                                    
                                                    %% Spectrum of the current collected from the emitter
                                                    j_e(n,:)=q...
                                                        *(f_ae(n,:)*eta_back_e_b(n,:)...
                                                        +f_be(n,:)*eta_forw_e_a(n,:)...
                                                        +f_ce(n,:)*eta_back_e_d(n,:)...
                                                        +f_de(n,:)*eta_forw_e_c(n,:));
                                                    
                                                    %% Spectrum of the current collected from the scr
                                                    j_scr(n,:)=q...
                                                        *(f_ascr(n,:)...
                                                        +f_bscr(n,:)...
                                                        +f_cscr(n,:)...
                                                        +f_dscr(n,:));
                                                    
                                                    %% Spectrum of the current collected from the base
                                                    j_b(n,:)=q...
                                                        *(f_ab(n,:)*eta_forw_b_a(n,:)...
                                                        +f_bb(n,:)*eta_back_b_b(n,:)...
                                                        +f_cb(n,:)*eta_forw_b_c(n,:)...
                                                        +f_db(n,:)*eta_back_b_d(n,:));
                                                    
                                                    EQE_e(n,:)=j_e(n,:)/(phi(n,:)*(h*c_0*1E9./lambda(n,:)))*1238/lambda(n,:);
                                                    EQE_scr(n,:)=j_scr(n,:)/(phi(n,:)*(h*c_0*1E9./lambda(n,:)))*1238/lambda(n,:);
                                                    EQE_b(n,:)=j_b(n,:)/(phi(n,:)*(h*c_0*1E9./lambda(n,:)))*1238/lambda(n,:);
                                                    EQE(n,:)=EQE_e(n,:)+EQE_scr(n,:)+EQE_b(n,:);
                                                    
                                                    IQE(n,:)=EQE(n,:)/A_si(n);
                                                    IQE_star(n,:)=EQE(n,:)/(1-R_tot(n,:));
                                                    n = n + 1;
                                                    
                                                    
                                                end
                                                
                                                if store_G==1
                                                    G_cell_tot=sum(G_cell)';
                                                    assignin('base', ['G_cell_' pol '_pol'] , G_cell_tot);
                                                end
                                                
                                                assignin('base', ['EQE_' pol '_pol'] , EQE_e+EQE_scr+EQE_b);
                                                assignin('base', ['IQE_' pol '_pol'] ,IQE);
                                                assignin('base', ['IQE_star_' pol '_pol'] ,IQE_star);
                                                assignin('base', ['R_tot_' pol '_pol'] ,R_tot);
                                                assignin('base', ['R_ext_' pol '_pol'] ,R_fe);
                                                assignin('base', ['R_esc_' pol '_pol'] ,R_esc);
                                                assignin('base', ['A_AR_' pol '_pol'] ,A_AR);
                                                assignin('base', ['A_si_' pol '_pol'] ,A_si);
                                                assignin('base', ['A_grid_' pol '_pol'] ,A_grid);
                                                assignin('base', ['A_back_tot_' pol '_pol'] ,A_back_tot);
                                                assignin('base', ['A_FCA_' pol '_pol'] ,A_FCA);
                                                assignin('base', ['r_b1_' pol '_pol'] ,r_b1);
                                                assignin('base', ['R_f1_' pol '_pol'] ,R_f1);
                                                
                                            end
                                            name_array(batch,:)={['theta: ' num2str(theta_0)] [', W: ' num2str(W)]...
                                                [', d_rel_b: ' num2str(d_rel_b)] [', d_rel_f: ' num2str(d_rel_f)]...
                                                [', r_db: ' num2str(r_db)] [', r_df: ' num2str(r_df)]...
                                                [', L_n: ' num2str(L_n)] [', L_p: ' num2str(L_p)]...
                                                [', S_e: ' num2str(S_e)] [', S_b: ' num2str(S_b)]};
                                            
                                            if store_G==1
                                            assignin('base', ['G_cell_' num2str(batch)], (G_cell_s_pol+G_cell_p_pol)/2);
                                            end
                                            
                                            assignin('base', ['EQE_' num2str(batch)], (EQE_s_pol+EQE_p_pol)/2);
                                            assignin('base', ['IQE_' num2str(batch)],(IQE_s_pol+IQE_p_pol)/2);
                                            assignin('base', ['IQE_star_' num2str(batch)],(IQE_star_s_pol+IQE_star_p_pol)/2);
                                            assignin('base', ['R_tot_' num2str(batch)],(R_tot_s_pol+R_tot_p_pol)/2);
                                            assignin('base', ['R_ext_' num2str(batch)],(R_ext_s_pol+R_ext_p_pol)/2);
                                            assignin('base', ['R_esc_' num2str(batch)],(R_esc_s_pol+R_esc_p_pol)/2);
                                            assignin('base', ['A_AR_' num2str(batch)],(A_AR_s_pol+A_AR_p_pol)/2);
                                            assignin('base', ['A_si_' num2str(batch)],(A_si_s_pol+A_si_p_pol)/2);
                                            assignin('base', ['A_grid_' num2str(batch)],(A_grid_s_pol+A_grid_p_pol)/2);
                                            assignin('base', ['A_back_tot_' num2str(batch)],(A_back_tot_s_pol+A_back_tot_p_pol)/2);
                                            assignin('base', ['A_FCA_' num2str(batch)],(A_FCA_s_pol+A_FCA_p_pol)/2);
                                            assignin('base', ['r_b1_' num2str(batch)],(r_b1_s_pol+r_b1_p_pol)/2);
                                            assignin('base', ['R_f1_' num2str(batch)],(R_f1_s_pol+R_f1_p_pol)/2);
                                            
                                            if batch==1
                                                theta_init=theta_0;
                                            elseif theta_0==theta_init
                                                batch_init=batch_init+1;
                                            end
                                            
                                            batch=batch+1
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
batch=batch-1;

