clear; clc; close all;

c = 3e8;                    
N_pa = 8;                    
P_subcarriers = 64;          

a_wg = 0.0055;                 
fc_cutoff = c / (2*a_wg);    

kappa0_design = 10;
f_design_pa = 30.5e9; 
coupling_efficiency_per_pa = 1;
kappa_L_target = asin(sqrt(coupling_efficiency_per_pa));
L_pa_length = kappa_L_target / kappa0_design;

delta_spacing = L_pa_length + c/f_design_pa/2; 

xu_user = 5; yu_user = 2; h_wg = 5;
xn_pa = xu_user + ((0:N_pa-1) - (N_pa-1)/2) * delta_spacing;

L_array = xn_pa(end) - xn_pa(1);

f_op_center_factors = [1.1, 1.5, 2.0];
f_op_center_values = f_op_center_factors * fc_cutoff;
BW_values_MHz = linspace(1, 2000, 20);
BW_values = BW_values_MHz * 1e6;
min_cp_overhead_percent = zeros(length(f_op_center_values), length(BW_values));

d_u_n = sqrt((xn_pa - xu_user).^2 + yu_user^2 + h_wg^2);
delta_tau_FS_array = (max(d_u_n) - min(d_u_n)) / c;


for i_f_op = 1:length(f_op_center_values)
    f_op_center = f_op_center_values(i_f_op);
    for i_bw = 1:length(BW_values)
        BW = BW_values(i_bw);
        T_sym_data_part = P_subcarriers / BW;
        f_edge1 = f_op_center - BW/2;
        f_edge2 = f_op_center + BW/2;
        
        if f_edge1 <= fc_cutoff
            min_cp_overhead_percent(i_f_op, i_bw) = NaN;
            continue;
        end
        
        vg_edge1 = c * sqrt(1 - (fc_cutoff/f_edge1)^2);
        vg_edge2 = c * sqrt(1 - (fc_cutoff/f_edge2)^2);
        delta_tau_g_max_waveguide = abs(1/vg_edge1 - 1/vg_edge2) * L_array;
        
        delta_tau_total = delta_tau_g_max_waveguide + delta_tau_FS_array;
        
        min_cp_overhead_percent(i_f_op, i_bw) = (delta_tau_total / T_sym_data_part) * 100;
    end
end

figure;
plot_markers_lines = {'-o', '-s', '-^'};
plot_colors = {'b', 'r', 'g'};
legend_entries = cell(1, length(f_op_center_values));
for i_f_op = 1:length(f_op_center_values)
    legend_entries{i_f_op} = sprintf('$f_{c} = %.1f f_{0}$', f_op_center_factors(i_f_op));
    plot(BW_values_MHz, min_cp_overhead_percent(i_f_op, :), plot_markers_lines{i_f_op}, ...
         'Color', plot_colors{i_f_op}, 'LineWidth', 1.5);
    hold on;
end
xlabel('System Bandwidth B (MHz)');
ylabel('Minimum Required CP Overhead (%)');
legend(legend_entries, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
grid off;
ax = gca;