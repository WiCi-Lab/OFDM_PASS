
clear; clc; close all;

c = 3e8; P = 64; f0 = 28e9; BW = 2e9;
delta_f = BW / P;
f_subcarriers = f0 - BW/2 + delta_f/2 : delta_f : f0 + BW/2;
Pt_dBm = 30; Pt = 10^((Pt_dBm-30)/10);
N0_dBm_Hz = -174; N0 = 10^((N0_dBm_Hz-30)/10); sigma2 = N0 * delta_f;
N = 8;
L_CP = P/4;

a = 0.01; fc = c / (2*a); h = 5; yu = 5; zu = 0;
kappa0 = 10; kappa1 = 2;
coupling_efficiency_per_pa = 1/2;
kappa_L_target = asin(sqrt(coupling_efficiency_per_pa));
L_pa = kappa_L_target / kappa0;
delta_x = L_pa + c/f0/2;
ne = 1.5;
delta_w = c/f0/ne/2;

loss_profiles = struct();
loss_profiles(1).name = 'No-Loss';
loss_profiles(1).coeffs.C1 = 0;
loss_profiles(1).coeffs.C2 = 0;
target_loss_low_Np = 0.1 / 8.686;
loss_profiles(2).name = '\alpha_{\rm{g}}=0.1dB/m';
loss_profiles(2).coeffs.C1 = 0.5 * target_loss_low_Np / sqrt(28e9);
loss_profiles(2).coeffs.C2 = 0.5 * target_loss_low_Np / 28e9;
target_loss_high_Np = 1.0 / 8.686;
loss_profiles(3).name = '\alpha_{\rm{g}}=1 dB/m';
loss_profiles(3).coeffs.C1 = 0.5 * target_loss_high_Np / sqrt(28e9);
loss_profiles(3).coeffs.C2 = 0.5 * target_loss_high_Np / 28e9;
matching_modes = {'perfectly_matched', 'mismatched'};

D_vector = linspace(10, 300, 20);
rate_pass = zeros(length(loss_profiles), length(matching_modes), length(D_vector));
rate_fixed_los = zeros(size(D_vector));
rate_fixed_nlos = zeros(size(D_vector)); 

params = struct('c',c,'P',P,'f0',f0,'fc',fc,'a',a,'h',h,'yu',yu,'zu',zu, ...
                'delta_x',delta_x,'kappa0',kappa0,'kappa1',kappa1,'L_pa',L_pa, ...
                'Pt',Pt,'P_p',Pt/P,'sigma2',sigma2,'f_subcarriers',f_subcarriers,'delta_f',delta_f,'BW',BW,'L_CP',L_CP);

standoff_dist = 5;
monte_carlo_runs = 200; 

for i_d = 1:length(D_vector)
    D = D_vector(i_d);
    user_pos = [D, yu, h];
    
    rate_fixed_los(i_d) = calculate_rate_fixed_ideal_miso(N, params, user_pos);
    
    rate_sum_nlos = 0;
    for mc_run = 1:monte_carlo_runs
        rate_sum_nlos = rate_sum_nlos + calculate_rate_fixed_nlos(N, params, user_pos);
    end
    rate_fixed_nlos(i_d) = rate_sum_nlos / monte_carlo_runs;
    
    disp(['Distance D = ', num2str(D), 'm complete.']);
end
disp('Benchmark calculations complete.');

for i_loss = 1:length(loss_profiles)
    alpha_coeffs = loss_profiles(i_loss).coeffs;
    sim_modes = matching_modes;
    if i_loss == 1
        sim_modes = {'perfectly_matched'};
    end
    for i_mode = 1:length(sim_modes)
        mode = sim_modes{i_mode};
        disp(['Simulating PASS: ', loss_profiles(i_loss).name, ', ', mode, '...']);
        for i_d = 1:length(D_vector)
            D = D_vector(i_d);
            pass_array_center_x = max(0, D - standoff_dist);
            user_pos = [D, yu, h];
            rate_pass(i_loss, i_mode, i_d) = calculate_rate_pass(N, params, alpha_coeffs, pass_array_center_x, user_pos, mode);
        end
    end

    if i_loss == 1
        sim_modes = {'perfectly_matched'};
        mode = sim_modes{i_mode};
        disp(['Simulating PASS: ', loss_profiles(i_loss).name, ', ', mode, '...']);
        for i_d = 1:length(D_vector)
            D = D_vector(i_d);
            pass_array_center_x = max(0, D - standoff_dist);
            user_pos = [D, yu, h];
            rate_pass(i_loss, i_mode, i_d) = calculate_ideal_rate_pass(N, params, alpha_coeffs, pass_array_center_x, user_pos, mode);
        end
    end

end
disp('All simulations complete.');

figure;
set(gcf, 'Position', [100, 100, 950, 700]);
hold on;
plot(D_vector, rate_fixed_los / 1e9, '-s', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Fixed-Location Antennas (LoS)');
plot(D_vector, rate_fixed_nlos / 1e9, '--s', 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', 'Fixed-Location Antennas (NLoS)');
plot(D_vector, squeeze(rate_pass(1, 1, :)) / 1e9, '-*', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'DisplayName', 'OFDM PASS (\alpha_{\rm{g}}=0, \Delta\beta = 0)');
colors = {'b', 'r'};
styles = {'-o', ':^'};
legend_labels = {'\Delta\beta = 0', '\Delta\beta \neq 0'};
for i_loss = 2:length(loss_profiles)
    for i_mode = 1:length(matching_modes)
        displayName = ['OFDM PASS (', loss_profiles(i_loss).name, ', ', legend_labels{i_mode}, ')'];
        plot(D_vector, squeeze(rate_pass(i_loss, i_mode, :)) / 1e9, ...
             styles{i_mode}, 'Color', colors{i_loss-1}, 'LineWidth', 2, ...
             'DisplayName', displayName);
    end
end
xlabel('User Longitudinal Distance from AP $x_u$ (m)', 'Interpreter', 'latex');
ylabel('Total Achievable Rate (Gbps)');
legend('show', 'Location', 'best','FontSize',13);
grid off;
box on
set(gca, 'FontSize', 13);
hold off;



function total_rate_bps = calculate_rate_fixed_nlos(N, params, user_pos)
    c=params.c; P=params.P; f_subcarriers=params.f_subcarriers;
    L_CP=params.L_CP; BW=params.BW;
    delta_f=params.delta_f; sigma2=params.sigma2; Pt=params.Pt;
    xu=user_pos(1); yu=user_pos(2); h=params.h;
    P_p = Pt / P;
    
    nlos_path_loss_exponent = 3.0;
    
    total_rate_bps_hz = 0;
    for p = 1:P
        fp = f_subcarriers(p);
        d_user = sqrt(xu^2 + yu^2 + h^2);
        
        path_loss_ref_db = 20*log10(4*pi*1*fp/c);
        path_loss_nlos_db = path_loss_ref_db + 10 * nlos_path_loss_exponent * log10(d_user);
        path_gain_nlos_linear = 10^(-path_loss_nlos_db / 10);
        
        h_fade = (randn(N, 1) + 1j * randn(N, 1)) / sqrt(2);
        h_nlos_total = sum(h_fade);
        
        channel_power_gain = path_gain_nlos_linear * abs(h_nlos_total)^2;
        SNR_p = (P_p / N) * channel_power_gain / sigma2;
        
        total_rate_bps_hz = total_rate_bps_hz + log2(1 + SNR_p);
    end
    total_rate_bps = total_rate_bps_hz* (BW / (P + L_CP));
end

function total_rate_bps = calculate_rate_fixed_ideal_miso(N, params, user_pos)
    c=params.c;P=params.P;f_subcarriers=params.f_subcarriers;delta_f=params.delta_f;
    
    L_CP=params.L_CP; BW=params.BW;
    sigma2=params.sigma2;Pt=params.Pt;xu=user_pos(1);yu=user_pos(2);h=params.h;
    P_p=Pt/P;total_rate_bps_hz=0;
    for p=1:P
        fp=f_subcarriers(p);d_user=sqrt(xu^2+yu^2+h^2);
        h_fs=(c/(4*pi*fp*d_user))*exp(-1j*2*pi*d_user*fp/c);
        channel_power_gain=N*abs(h_fs)^2;
        SNR_p=P_p*channel_power_gain/sigma2;
        total_rate_bps_hz=total_rate_bps_hz+log2(1+SNR_p);
    end
    total_rate_bps=total_rate_bps_hz* (BW / (P + L_CP));
end

function total_rate_bps = calculate_rate_pass(N, params, alpha_coeffs, array_center_x, user_pos, mode)
    c=params.c;P=params.P;f0_center=params.f0;fc=params.fc;a=params.a;h=params.h;
    delta_x=params.delta_x;L_pa=params.L_pa;P_p=params.P_p;sigma2=params.sigma2;
    L_CP=params.L_CP; BW=params.BW;
    f_subcarriers=params.f_subcarriers;delta_f=params.delta_f;kappa0=params.kappa0;
    kappa1=params.kappa1;C1=alpha_coeffs.C1;C2=alpha_coeffs.C2;
    xu=user_pos(1);yu=user_pos(2);zu=0;
    get_kappa=@(fp)kappa0+kappa1*((fp-f0_center)/f0_center);
    get_alpha_g=@(fp)C1*sqrt(fp)+C2*fp;
    get_beta_g=@(fp)sqrt((2*pi*fp/c).^2-(pi/a)^2);
    beta_g_at_f0=get_beta_g(f0_center);
    np_center_matched=beta_g_at_f0/(2*pi*f0_center/c);
    xn=array_center_x+((0:N-1)-(N-1)/2)*delta_x;
    total_rate_bps_hz=0;

    
    for p=1:P
        fp=f_subcarriers(p);if fp<=fc,continue;end
        alpha_g_f=get_alpha_g(fp);beta_g_f=get_beta_g(fp);kappa_f=get_kappa(fp);
        initial_propagation_dist=xn(1);if initial_propagation_dist<0,initial_propagation_dist=0;end
        T_initial=exp(-(alpha_g_f+1j*beta_g_f)*initial_propagation_dist);
        H_p_sum=0;E_wg_rem=T_initial;
        for n=1:N
            E_in_n=E_wg_rem;
            if strcmp(mode,'perfectly_matched'),delta_beta_f=0;else,np_mismatched=1.5;beta_p_f=(2*pi*fp/c)*np_mismatched;delta_beta_f=beta_g_f-beta_p_f;end
            S_p_n=sqrt(kappa_f^2+(delta_beta_f/2)^2);
            alpha_coup_n=-1j*(kappa_f/S_p_n)*sin(S_p_n*L_pa)*exp(1j*delta_beta_f*L_pa/2);
            E_coup_n=E_in_n*alpha_coup_n;
            d_un=sqrt((xn(n)-xu)^2+yu^2+(h-zu)^2);
            h_n_p=(c/(4*pi*fp*d_un))*exp(-1j*2*pi*d_un*fp/c);
            H_p_sum=H_p_sum+abs(E_coup_n*h_n_p);
            alpha_wg_n=(cos(S_p_n*L_pa)+1j*(delta_beta_f/(2*S_p_n))*sin(S_p_n*L_pa))*exp(-1j*delta_beta_f*L_pa/2);
            E_after_pa_n=E_in_n*alpha_wg_n;
            if n<N,T_segment=exp(-(alpha_g_f+1j*beta_g_f)*(xn(n+1)-xn(n)));E_wg_rem=E_after_pa_n*T_segment;end
        end
        SNR_p=P_p*abs(H_p_sum)^2/sigma2;
        total_rate_bps_hz=total_rate_bps_hz+log2(1+SNR_p);
    end
    total_rate_bps=total_rate_bps_hz* (BW / (P + L_CP));
end

function total_rate_bps = calculate_ideal_rate_pass(N, params, alpha_coeffs, array_center_x, user_pos, mode)
    c=params.c;P=params.P;f0_center=params.f0;fc=params.fc;a=params.a;h=params.h;
    delta_x=params.delta_x;L_pa=params.L_pa;P_p=params.P_p;sigma2=params.sigma2;
    L_CP=params.L_CP; BW=params.BW;
    f_subcarriers=params.f_subcarriers;delta_f=params.delta_f;kappa0=params.kappa0;
    kappa1=params.kappa1;C1=alpha_coeffs.C1;C2=alpha_coeffs.C2;
    xu=user_pos(1);yu=user_pos(2);zu=0;
    get_kappa=@(fp)kappa0+kappa1*((fp-f0_center)/f0_center);
    get_alpha_g=@(fp)C1*sqrt(fp)+C2*fp;
    get_beta_g=@(fp)sqrt((2*pi*fp/c).^2-(pi/a)^2);
    beta_g_at_f0=get_beta_g(f0_center);
    np_center_matched=beta_g_at_f0/(2*pi*f0_center/c);
    xn=array_center_x+((0:N-1)-(N-1)/2)*delta_x;
    total_rate_bps_hz=0;

    total_radiated_power_fraction = 1 - (1 - 1/N)^N;
    power_compensation_factor = 1 / total_radiated_power_fraction;
    
    for p=1:P
        fp=f_subcarriers(p);if fp<=fc,continue;end
        alpha_g_f=get_alpha_g(fp);beta_g_f=get_beta_g(fp);kappa_f=get_kappa(fp);
        initial_propagation_dist=xn(1);if initial_propagation_dist<0,initial_propagation_dist=0;end
        T_initial=exp(-(alpha_g_f+1j*beta_g_f)*initial_propagation_dist);
        H_p_sum=0;E_wg_rem=T_initial;
        for n=1:N
            E_in_n=E_wg_rem;
            if strcmp(mode,'perfectly_matched'),delta_beta_f=0;else,np_mismatched=1.5;beta_p_f=(2*pi*fp/c)*np_mismatched;delta_beta_f=beta_g_f-beta_p_f;end
            S_p_n=sqrt(kappa_f^2+(delta_beta_f/2)^2);
            alpha_coup_n=-1j*(kappa_f/S_p_n)*sin(S_p_n*L_pa)*exp(1j*delta_beta_f*L_pa/2);
            E_coup_n=E_in_n*alpha_coup_n;
            d_un=sqrt((xn(n)-xu)^2+yu^2+(h-zu)^2);
            h_n_p=(c/(4*pi*fp*d_un))*exp(-1j*2*pi*d_un*fp/c);
            H_p_sum=H_p_sum+abs(E_coup_n*h_n_p);
            alpha_wg_n=(cos(S_p_n*L_pa)+1j*(delta_beta_f/(2*S_p_n))*sin(S_p_n*L_pa))*exp(-1j*delta_beta_f*L_pa/2);
            E_after_pa_n=E_in_n*alpha_wg_n;
            if n<N,T_segment=exp(-(alpha_g_f+1j*beta_g_f)*(xn(n+1)-xn(n)));E_wg_rem=E_after_pa_n*T_segment;end
        end
        SNR_p=P_p*abs(H_p_sum)^2/sigma2*power_compensation_factor;
        total_rate_bps_hz=total_rate_bps_hz+log2(1+SNR_p);
    end
    total_rate_bps=total_rate_bps_hz* (BW / (P + L_CP));
end