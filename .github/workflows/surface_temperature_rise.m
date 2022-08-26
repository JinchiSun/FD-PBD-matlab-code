function Delta_Ts = surface_temperature_rise(ff)

% from top to bottome: 1: air: 2: Al; 3: bulk material

w0 = 8.3e-6; % pump and probe beam radii
A0 = 0.282e-3; % absorbed power of modulated pump beam

% 1: air
sigma_1 = 0.028;
capac_1 = 1192;
Dif_1 = sigma_1/capac_1;

% 2: Al
L_2 = 85e-9;
sigma_2 = 165;
capac_2 = 2.42e6;
Dif_2 = sigma_2/capac_2;
G_int = 40e6;

% 3: (100) SrTiO3
sigma_3_z = 11;
sigma_3_r_z = 1; % sigma_3_r/sigma_3_z
capac_3 = 2.74e6;
Dif_3 = sigma_3_z/capac_3;
%%
n_p = 3001;
up_p = 8/w0;
d_p = up_p/n_p;
pp = d_p:d_p:up_p;

avg_thetas_p_omega = zeros(n_p,length(ff));
Delta_Ts = ones(length(ff),1);

for i_fr = 1:length(ff)
omega = 2*pi*ff(i_fr);
    qn2_1 = (1i)*omega/Dif_1; % air
    qn2_2 = (1i)*omega/Dif_2; % Al
    qn2_3 = (1i)*omega/Dif_3; % bulk material
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
        avg_thetas_p_omega(i_p,i_fr) = 1/(2*pi)*theta_s*exp(-w0^2*p^2/8)*p;
    end
    Delta_Ts(i_fr,1) = simpson_inte(avg_thetas_p_omega(:,i_fr),d_p);
end

end

function I = simpson_inte(array,pace)
steps = length(array);
edge_sum = sum(array(3:2:steps-2));
mid_sum = sum(array(2:2:steps-1));
I = (1/6)*(2*pace)*(array(1) + 2*edge_sum + 4*mid_sum + array(steps));
end
