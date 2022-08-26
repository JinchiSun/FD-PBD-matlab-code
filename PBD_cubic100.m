function PBD_ani_r0_phi = PBD_cubic100(ff)

% from top to bottome: 1: air: 2: Al; 3: bulk material
% PBD_ani_r0 = PBD calculated using the full model
% for PBD due to surface deformation of Al coated bulk material sample
% the bulk material is cubic (100)

w0 = 8.3e-6; % pump and probe beam radii
r0 = 9.8e-6; % beam offset distance
phi = 0/180*pi; % beam offset direction vector relative to lattice
A0 = 0.282e-3; % absorbed power of modulated pump beam

% 1: air
sigma_1 = 0.028;
capac_1 = 1192;
Dif_1 = sigma_1/capac_1;

% 2: Al (with (111) texture)
L_2 = 85e-9;
sigma_2 = 165;
capac_2 = 2.42e6;
Dif_2 = sigma_2/capac_2;
G_int = 40e6;

C11_3 = 107.4e9;
C12_3 = 60.5e9;
C44_3 = 28.3e9;
rho_2 = 2.70e3;
alphaT_3 = 23.1e-6;

C11_2_p = (C11_3 + C12_3 + 2*C44_3)/2;
C33_2_p = (C11_3 + 2*C12_3 + 4*C44_3)/3;
C44_2_p = (C11_3 - C12_3 + C44_3)/3;
C12_2_p = (C11_3 + 5*C12_3 - 2*C44_3)/6;
C13_2_p = (C11_3 + 2*C12_3 - 4*C44_3)/3;
C46_2_p = 0;
C22_2_p = C11_2_p;
C23_2_p = C13_2_p;
C55_2_p = C44_2_p;
C66_2_p = (C11_2_p - C12_2_p)/2;
beta_3 = (C11_3 + 2*C12_3)*alphaT_3;
betax_2 = beta_3;
betay_2 = beta_3;
betaz_2 = beta_3;

C22C11_2 = C22_2_p/C11_2_p;
C33C11_2 = C33_2_p/C11_2_p;
C12C11_2 = C12_2_p/C11_2_p;
C13C11_2 = C13_2_p/C11_2_p;
C23C11_2 = C23_2_p/C11_2_p;
C44C11_2 = C44_2_p/C11_2_p;
C55C11_2 = C55_2_p/C11_2_p;
C66C11_2 = C66_2_p/C11_2_p;
C46C11_2 = C46_2_p/C11_2_p;
sqrtC11rho_2 = sqrt(C11_2_p/rho_2);
betaxC11_2 = betax_2/C11_2_p;
betayC11_2 = betay_2/C11_2_p;
betazC11_2 = betaz_2/C11_2_p;


% 3: SrTiO3 (100)
sigma_3_z = 11;
sigma_3_r_z = 1; % sigma_3_r/sigma_3_z
capac_3 = 2.74e6;
Dif_3 = sigma_3_z/capac_3;

C11_3 = 317.2e9;
C12_3 = 102.5e9;
C44_3 = 123.5e9;
alphaT_3 = 10.4e-6;
rho_3 = 5.175e3;

C12C11_3 = C12_3/C11_3;
C44C11_3 = C44_3/C11_3;
sqrtC11rho_3 = sqrt((1+(1e-6)*(1i))*C11_3/rho_3);
beta_3 = (C11_3 + 2*C12_3)*alphaT_3;
betaC11_3 = beta_3/C11_3;


%%
n_p = 3001;
up_p = 8/w0;
d_p = up_p/n_p;
pp = d_p:d_p:up_p;

n_psi = 45;
up_psi = pi/4;
d_psi = up_psi/n_psi;
ppsi = 0:d_psi:(up_psi - d_psi);

C_probe = 1;
Z_p_psi_omega = zeros(n_p,n_psi,length(ff));
PBD_ani_r0_phi = ones(length(ff),1);

A_2 = zeros(6,6);
B_2 = zeros(6,6);
D_2 = zeros(6,1);
A_2(1,4) = 1;
A_2(2,5) = 1;
A_2(3,6) = 1;
A_2(4,1) = C55C11_2;
A_2(5,2) = C44C11_2;
A_2(6,3) = C33C11_2;
B_2(4,4) = 1;
B_2(5,5) = 1;
B_2(6,6) = 1;
D_2(6) = betazC11_2;

A_3 = zeros(6,6);
B_3 = zeros(6,6);
D_3 = zeros(6,1);
A_3(1,4) = 1;
A_3(2,5) = 1;
A_3(3,6) = 1;
A_3(4,1) = C44C11_3;
A_3(5,2) = C44C11_3;
A_3(6,3) = 1;
B_3(4,4) = 1;
B_3(5,5) = 1;
B_3(6,6) = 1;
D_3(6) = betaC11_3;

for i_fr = 1:length(ff)
    omega = 2*pi*ff(i_fr);
    qn2_1 = (1i)*omega/Dif_1;
    qn2_2 = (1i)*omega/Dif_2;
    qn2_3 = (1i)*omega/Dif_3;
    I_p2 = zeros(n_p,1);
    for i_p = 1:n_p
        p = pp(i_p);
        flx = A0*exp(-w0^2*p^2/8);
        zeta_1 = sqrt(qn2_1 + p^2);
        zeta_2 = sqrt(qn2_2 + p^2);
        zeta_3 = sqrt(qn2_3 + sigma_3_r_z*p^2);
        zeta_2L_2 = zeta_2*L_2;
        sigma_1zeta_1 = sigma_1*zeta_1;
        sigma_2zeta_2 = sigma_2*zeta_2;
        sigma_3zeta_3 = sigma_3_z*zeta_3;
        G_d =  (sigma_3zeta_3*sinh(zeta_2L_2) + sigma_2zeta_2*cosh(zeta_2L_2) + sigma_2zeta_2*sigma_3zeta_3/G_int*cosh(zeta_2L_2))/sigma_2zeta_2;
        G_d = G_d/(sigma_3zeta_3*cosh(zeta_2L_2) + sigma_2zeta_2*sinh(zeta_2L_2) + sigma_2zeta_2*sigma_3zeta_3/G_int*sinh(zeta_2L_2));
        G_u = 1/sigma_1zeta_1;
        G = 1/(1/G_u + 1/G_d);
        theta_s = flx*G;
        theta_bs = cosh(zeta_2L_2)*theta_s + sigma_2zeta_2/G_int*sinh(zeta_2L_2)*theta_s - sinh(zeta_2L_2)*flx/sigma_2zeta_2 - cosh(zeta_2L_2)*flx/G_int;
        a_minus = (sigma_3zeta_3/G_int*theta_bs + theta_bs - theta_s*exp(-zeta_2L_2))/(exp(zeta_2L_2) - exp(-zeta_2L_2));
        a_plus = theta_s - a_minus;
        for i_psi = 1:n_psi
            psi = ppsi(i_psi);
            k = p*cos(psi);
            xi = p*sin(psi);

            A_2(1,1) = -C46C11_2*(1i)*k;
            A_2(2,2) = C46C11_2*(1i)*k;
            A_2(1,2) = C46C11_2*(1i)*xi;
            A_2(2,1) = C46C11_2*(1i)*xi;
            A_2(1,3) = C13C11_2*(1i)*k;
            A_2(2,3) = C23C11_2*(1i)*xi;
            B_2(1,1) = k^2 + C66C11_2*xi^2 - omega^2/sqrtC11rho_2^2;
            B_2(2,2) = C22C11_2*xi^2 + C66C11_2*k^2 - omega^2/sqrtC11rho_2^2;
            B_2(1,2) = (C12C11_2 + C66C11_2)*k*xi;
            B_2(2,1) = (C12C11_2 + C66C11_2)*k*xi;
            B_2(3,3) = -omega^2/sqrtC11rho_2^2;
            B_2(1,3) = C46C11_2*(xi^2 - k^2);
            B_2(2,3) = 2*C46C11_2*k*xi;
            B_2(3,4) = -(1i)*k;
            B_2(3,5) = -(1i)*xi;
            B_2(4,1) = C46C11_2*(1i)*k;
            B_2(4,2) = -C46C11_2*(1i)*xi;
            B_2(4,3) = -C55C11_2*(1i)*k;
            B_2(5,1) = -C46C11_2*(1i)*xi;
            B_2(5,2) = -C46C11_2*(1i)*k;
            B_2(5,3) = -C44C11_2*(1i)*xi;
            B_2(6,1) = -C13C11_2*(1i)*k;
            B_2(6,2) = -C23C11_2*(1i)*xi;
            D_2(1) = betaxC11_2*(1i)*k;
            D_2(2) = betayC11_2*(1i)*xi;
            inv_A_2 = inv(A_2);
            N_2 = inv_A_2*D_2;
            [Q_2, R_2] = eig(B_2,A_2);
            LAMBDA_2 = zeros(6,1);
            for i = 1:6
                LAMBDA_2(i) = R_2(i,i);
            end
            condi_2 = rcond(Q_2);
            if condi_2 < 2e-16
                pause;
            end
            inv_Q_2 = inv(Q_2);
            U_2 = inv_Q_2*N_2;

            A_3(1,3) = C12C11_3*(1i)*k;
            A_3(2,3) = C12C11_3*(1i)*xi;
            B_3(1,1) = k^2 + C44C11_3*xi^2 - omega^2/sqrtC11rho_3^2;
            B_3(2,2) = xi^2 + C44C11_3*k^2 - omega^2/sqrtC11rho_3^2;
            B_3(1,2) = (C12C11_3 + C44C11_3)*k*xi;
            B_3(2,1) = (C12C11_3 + C44C11_3)*k*xi;
            B_3(3,3) = -omega^2/sqrtC11rho_3^2;
            B_3(3,4) = -(1i)*k;
            B_3(3,5) = -(1i)*xi;
            B_3(4,3) = -C44C11_3*(1i)*k;
            B_3(5,3) = -C44C11_3*(1i)*xi;
            B_3(6,1) = -C12C11_3*(1i)*k;
            B_3(6,2) = -C12C11_3*(1i)*xi;
            D_3(1) = betaC11_3*(1i)*k;
            D_3(2) = betaC11_3*(1i)*xi;
            inv_A_3 = inv(A_3);
            N_3 = inv_A_3*D_3;
            [Q_raw_3, R_raw_3] = eig(B_3,A_3);
            Q_3 = zeros(6,6);
            R_3 = zeros(6,6);
            count_3 = 0;
            for i = 1:6
                if real(R_raw_3(i,i)) < 0
                    count_3 = count_3 + 1;
                    R_3(count_3,count_3) = R_raw_3(i,i);
                    Q_3(:,count_3) = Q_raw_3(:,i);
                else
                    R_3(i-count_3+3,i-count_3+3) = R_raw_3(i,i);
                    Q_3(:,i-count_3+3) = Q_raw_3(:,i);
                end
            end
            LAMBDA_3 = zeros(6,1);
            for i = 1:6
                LAMBDA_3(i) = R_3(i,i);
            end
            condi_3 = rcond(Q_3);
            if condi_3 < 2e-16
                pause;
            end
            inv_Q_3 = inv(Q_3);
            U_3 = inv_Q_3*N_3;

            BCM = zeros(9,9);
            BCC = zeros(9,1);
            for i_cl = 1:6
                BCM(1:3,i_cl) = Q_2(4:6,i_cl);
                BCM(4:9,i_cl) = Q_2(1:6,i_cl)*exp(LAMBDA_2(i_cl)*L_2);
            end
            for i_cl = 7:9
                BCM(4:6,i_cl) = -Q_3(1:3,i_cl-6)*exp(LAMBDA_3(i_cl-6)*L_2);
                BCM(7:9,i_cl) = -C11_3/C11_2_p*Q_3(4:6,i_cl-6)*exp(LAMBDA_3(i_cl-6)*L_2);
            end
            for i_rw = 1:3
                sum = 0;
                for jj = 1:6
                    temp = Q_2(i_rw+3,jj)*U_2(jj)*(a_minus/(zeta_2-LAMBDA_2(jj)) + a_plus/(-zeta_2-LAMBDA_2(jj)));
                    sum = sum + temp;
                end
                BCC(i_rw) = -sum;
            end
            for i_rw = 4:6
                sum1 = 0;
                sum2 = 0;
                for jj = 1:6
                    temp1 = Q_2(i_rw-3,jj)*U_2(jj)*(a_minus*exp(zeta_2L_2)/(zeta_2-LAMBDA_2(jj)) + a_plus*exp(-zeta_2L_2)/(-zeta_2-LAMBDA_2(jj)));
                    temp2 = Q_3(i_rw-3,jj)*U_3(jj)*(theta_bs/(-zeta_3-LAMBDA_3(jj)));
                    sum1 = sum1 + temp1;
                    sum2 = sum2 + temp2;
                end
                BCC(i_rw) = -sum1 + sum2;
            end
            for i_rw = 7:9
                sum1 = 0;
                sum2 = 0;
                for jj = 1:6
                    temp1 = Q_2(i_rw-3,jj)*U_2(jj)*(a_minus*exp(zeta_2L_2)/(zeta_2-LAMBDA_2(jj)) + a_plus*exp(-zeta_2L_2)/(-zeta_2-LAMBDA_2(jj)));
                    temp2 = Q_3(i_rw-3,jj)*U_3(jj)*(theta_bs/(-zeta_3-LAMBDA_3(jj)));
                    sum1 = sum1 + temp1;
                    sum2 = sum2 + temp2;
                end
                BCC(i_rw) = -sum1 + (C11_3/C11_2_p)*sum2;
            end
            J = BCM\BCC;
            w_s_H = Q_2(3,1)*J(1) + Q_2(3,2)*J(2) + Q_2(3,3)*J(3) + Q_2(3,4)*J(4) + Q_2(3,5)*J(5) + Q_2(3,6)*J(6);
            sssum = 0;
            for jj = 1:6
                temp = Q_2(3,jj)*U_2(jj)*(a_minus/(zeta_2-LAMBDA_2(jj)) - a_plus/(zeta_2+LAMBDA_2(jj)));
                sssum = sssum + temp;
            end
            w_s_P = sssum;
            Z_p_psi_omega(i_p,i_psi,i_fr) = -(w_s_H + w_s_P);
        end
    end
end

for i_fr = 1:length(ff)
    I_p2 = zeros(n_p,1);
    for i_p = 1:n_p
        p = pp(i_p);
        pr0 = p*r0;
        Zg = zeros(n_psi,1);
        for i_psi = 1:n_psi
            psi = ppsi(i_psi);
            I1 = I_part(phi,psi,pr0);
            I2 = I_part(phi,pi/2-psi,pr0);
            I3 = I_part(phi,pi/2+psi,pr0);
            I4 = I_part(phi,pi-psi,pr0);
            I5 = I_part(phi,pi+psi,pr0);
            I6 = I_part(phi,3*pi/2-psi,pr0);
            I7 = I_part(phi,3*pi/2+psi,pr0);
            I8 = I_part(phi,2*pi-psi,pr0);
            Zg(i_psi) = Z_p_psi_omega(i_p,i_psi,i_fr)*(I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8);
        end
        I_p2(i_p) = 1/(2*pi)*simpson_inte(Zg,d_psi)*exp(-w0^2*p^2/8)*p^2;
    end
    PBD_ani_r0_phi(i_fr,1) = C_probe/pi*simpson_inte(I_p2,d_p);
end

end

function Fpsiphi = I_part(phi,psi,pr0)
dif = psi-phi;
Fpsiphi = (1i)*cos(dif)*exp((1i)*pr0*cos(dif));
end

function I = simpson_inte(array,pace)
steps = length(array);
edge_sum = sum(array(3:2:steps-2));
mid_sum = sum(array(2:2:steps-1));
I = (1/6)*(2*pace)*(array(1) + 2*edge_sum + 4*mid_sum + array(steps));
end