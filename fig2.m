clear; clc; close all;

c = 3e8;                        % Speed of light (m/s)
N = 8;                          % Number of PAs per waveguide
P = 64;                         % Number of OFDM subcarriers
f0 = 28e9;                      % Center frequency (Hz)
BW = 2e9;                       % System bandwidth (Hz)
delta_f = BW / P;               % Subcarrier spacing (Hz)
f_subcarriers = f0 - BW/2 + delta_f/2 : delta_f : f0 + BW/2 - delta_f/2; %subcarrier frequencies

Pt_dBm = 30;                    % Total transmit power (dBm)
Pt = 10^((Pt_dBm - 30)/10);      % Total transmit power (W)
P_p = Pt / P;                   % Power per subcarrier
N0_dBm_Hz = -174;               % Noise power spectral density (dBm/Hz)
N0 = 10^((N0_dBm_Hz - 30)/10);   % Noise power spectral density (W/Hz)
sigma2 = N0 * delta_f;          % Noise power per subcarrier

a = 0.0055;                     % Waveguide width (m)
fc = c / (2*a);                 % Cutoff frequency (Hz) - TE10 mode
h = 5;                          % Waveguide height (m)

ap_x = 0; % AP (signal source) located at x=0
xu = 5; yu = 2; zu = 0;          % User coordinates (m)

kappa0 = 10;                    % Coupling coefficient at center frequency f0 (1/m)
kappa1 = 5;                     % Linear variation slope for kappa vs. frequency


params_alpha.C1 = 1e-6; % conductor loss
params_alpha.C2 = 6e-12; % dielectric loss

np_val = 1.5;                   % Effective refractive index of PA (assumed constant for simplicity)


L_pa = asin(sqrt(1/N)) / kappa0; % PA coupling length

delta_x = c/f0/2; % Spacing between PAs (half-wavelength at center frequency)
xn_approx = xu + ((0:N-1) - (N-1)/2) * delta_x;

Rp_freq_dependent = zeros(1, P);
Rp_freq_independent = zeros(1, P);
Rp_delta_beta_zero = zeros(1, P);

params = struct('c', c, 'N', N, 'f0', f0, 'fc', fc, 'a', a, 'h', h, ...
                'ap_x', ap_x, 'xu', xu, 'yu', yu, 'zu', zu, ...
                'kappa0', kappa0, 'kappa1', kappa1, 'L_pa', L_pa, ...
                'np_val', np_val, 'params_alpha', params_alpha, ...
                'xn_approx', xn_approx, 'P_p', P_p, 'sigma2', sigma2);

for p = 1:P
    fp = f_subcarriers(p);
    % Case 1: Practical model where parameters are frequency dependent
    Rp_freq_dependent(p) = calc_channel_revised(fp, params, 'dependent');
    % Case 2: Simplified model where parameters are fixed at center frequency
    Rp_freq_independent(p) = calc_channel_revised(fp, params, 'independent');
    % Case 3: Ideal model with perfect phase matching (delta_beta = 0)
    Rp_delta_beta_zero(p) = calc_channel_revised(fp, params, 'delta_beta_zero');
end

figure;
h1 = plot(f_subcarriers/1e9, Rp_delta_beta_zero, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Ideal Freq. Dependent ($\Delta\beta = 0$)');
hold on;
h2 = plot(f_subcarriers/1e9, Rp_freq_dependent, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Practical Freq. Dependent');
h3 = plot(f_subcarriers/1e9, Rp_freq_independent, 'g-^', 'LineWidth', 1.5, 'DisplayName', 'Freq. Independent');

fc_ghz = fc / 1e9;
current_ylim = ylim;
line([fc_ghz fc_ghz], current_ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
label_text = {sprintf('Cutoff frequency:'); sprintf('f_0 = %.1f GHz', fc_ghz)};
text(fc_ghz*1.005, mean(current_ylim), label_text, 'HorizontalAlignment', 'left', 'FontSize', 10);

xlabel('Subcarrier Frequency (GHz)');
ylabel('Subcarrier Rate (bps/Hz)');
legend([h1, h2, h3], 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
box on;
hold off;


%% --- Channel Calculation Function ---
function Rp = calc_channel_revised(fp, params, mode)
    c = params.c; N = params.N; f0_center = params.f0; fc = params.fc; a = params.a;
    h = params.h; ap_x = params.ap_x; xu = params.xu; yu = params.yu; zu = params.zu;
    L_pa = params.L_pa; np_val = params.np_val; xn_approx = params.xn_approx;
    P_p = params.P_p; sigma2 = params.sigma2;
    kappa0 = params.kappa0; kappa1 = params.kappa1;
    params_alpha = params.params_alpha;

    get_kappa = @(freq) kappa0 + kappa1 * ((freq - f0_center) / f0_center);
    get_alpha_g = @(freq) params_alpha.C1 * sqrt(freq) + params_alpha.C2 * freq;
    get_beta_g = @(freq) sqrt((2*pi*freq/c).^2 - (pi/a)^2);

    freq_calc = fp;
    if strcmp(mode, 'independent')
        % All physical parameters are calculated at the center frequency
        freq_calc = f0_center;
    end
    
    % If the signal frequency is below the cutoff frequency, it cannot propagate
    if fp <= fc
        H_p_sum = 0;
    else
        % Calculate frequency-dependent parameters for this subcarrier
        beta_g_f = get_beta_g(freq_calc);      % Waveguide propagation constant
        alpha_g_f = get_alpha_g(freq_calc);     % Waveguide attenuation constant
        gamma_f = alpha_g_f + 1j * beta_g_f;    % Complex propagation constant
        kappa_f = get_kappa(freq_calc);         % Coupling coefficient
        beta_p_f = (2*pi*freq_calc/c) * np_val; % PA phase constant
        delta_beta_f = beta_g_f - beta_p_f;     % Phase mismatch

        % Force phase mismatch to zero for ideal cases
        if strcmp(mode, 'delta_beta_zero') || strcmp(mode, 'independent')
            delta_beta_f = 0; 
        end

        H_p_sum = 0;
        
        % Calculate the initial propagation from the source (ap_x=0) to the first PA
        initial_prop_dist = xn_approx(1) - ap_x;
        T_initial = exp(-gamma_f * initial_prop_dist);
        
        % Complex amplitude at the entrance of the first PA coupling region
        E_wg_rem = T_initial; 

        for n = 1:N
            % Complex amplitude of the signal at the entrance of the n-th PA
            E_in_n = E_wg_rem;
            
            % Coupled-mode theory parameter S_p
            S_p_n = sqrt(kappa_f^2 + (delta_beta_f / 2)^2);
            
            % Local coupling factor (describes power transfer from waveguide to PA)
            alpha_coup_n = -1j * (kappa_f / S_p_n) * sin(S_p_n * L_pa) * exp(1j * delta_beta_f * L_pa / 2);
            
            % Free-space channel from PA n to user u (using actual subcarrier frequency fp)
            d_un = sqrt((xn_approx(n) - xu)^2 + yu^2 + (h - zu)^2);
            h_n_p = (c / (4*pi*fp*d_un)) * exp(-1j * 2*pi * d_un * fp / c);
            
            % Accumulate the total effective channel contribution from PA n
            H_p_sum = H_p_sum + E_in_n * alpha_coup_n * h_n_p;
            
            % Signal remaining in the waveguide after coupling through PA n (through-factor)
            alpha_wg_n = (cos(S_p_n * L_pa) + 1j * (delta_beta_f / (2 * S_p_n)) * sin(S_p_n * L_pa)) * exp(-1j * delta_beta_f * L_pa / 2);
            E_after_coupling = E_in_n * alpha_wg_n;
            
            % The remaining signal propagates from PA n to PA n+1
            if n < N
                inter_pa_dist = xn_approx(n+1) - xn_approx(n);
                T_segment = exp(-gamma_f * inter_pa_dist);
                E_wg_rem = E_after_coupling * T_segment;
            end
        end
    end
    
    SNR_p = P_p * abs(H_p_sum)^2 / sigma2;
    Rp = log2(1 + SNR_p);
end
