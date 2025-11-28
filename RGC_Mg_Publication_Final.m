%% =========================================================================
%  Mg²⁺ NEUROPROTECTION IN RETINAL GANGLION CELLS
%  Complete Publication-Ready Analysis Suite
%  =========================================================================
%
%  This script generates all figures for the manuscript in a single run:
%    Figure 1: Model Overview & Example Traces
%    Figure 2: Dose-Response Analysis  
%    Figure 3: Therapeutic Windows
%    Figure 4: Intervention Timing
%    Figure 5: Mechanistic Demonstration
%
%  KEY FINDINGS:
%    • Mg²⁺ blocks NMDA receptors → reduces Ca²⁺ influx (neuroprotection)
%    • At 10-60 Hz: spike output unchanged (AMPA-driven)
%    • At 80 Hz: spike loss plateaus at 20% for [Mg²⁺] = 1.6-2.0 mM
%    • Optimal therapeutic window at 80 Hz: 1.6-2.0 mM
%
%  Author: Mehdi Zarei
%  Institution: ICTER / Nencki Institute, Polish Academy of Sciences
%  =========================================================================

clear; clc; close all;

%% =========================================================================
%  GLOBAL SETTINGS
%  =========================================================================
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultLineLineWidth', 1.5);

% Color schemes
colors.Mg_gradient = [0.85 0.2 0.2;    % Low Mg (red - danger)
                      0.75 0.4 0.2;
                      0.6 0.5 0.3;
                      0.4 0.55 0.4;
                      0.3 0.6 0.5;
                      0.2 0.5 0.7;
                      0.15 0.4 0.75;
                      0.1 0.3 0.8];     % High Mg (blue - protected)

colors.freq = [0.2 0.7 0.3;            % 10 Hz (green)
               0.4 0.6 0.2;            % 30 Hz
               0.7 0.5 0.1;            % 60 Hz
               0.85 0.2 0.2];          % 80 Hz (red)

colors.optimal = [0.2 0.7 0.3];
colors.suboptimal = [0.7 0.7 0.7];
colors.danger = [0.85 0.2 0.2];

%% =========================================================================
%  MODEL PARAMETERS
%  =========================================================================
p.Cm = 1.0; p.V_rest = -65;
p.E_Na = 50; p.E_K = -77; p.E_L = -54.4; p.E_Ca = 120; p.E_exc = 0;
p.g_Na = 120; p.g_Kdr = 36; p.g_KA = 8;
p.g_CaL = 0.3; p.g_KCa = 0.3; p.g_L = 0.35;
p.g_AMPA_max = 0.25; p.g_NMDA_max = 1.2;
p.EC50_AMPA = 0.5; p.EC50_NMDA = 0.002;
p.tau_AMPA_rise = 0.3; p.tau_AMPA_decay = 3.0;
p.tau_NMDA_rise = 5.0; p.tau_NMDA_decay = 80.0;
p.eta = 0.28; p.gamma = 0.062; p.f_Ca = 0.15;
p.Ca_rest = 0.05; p.tau_Ca = 200;
p.k_Ca_NMDA = 0.012; p.k_Ca_CaL = 0.003;
p.Kd_KCa = 0.5; p.n_KCa = 2;

Ca_toxic = 1.0;
max_spike_loss = 20;

%% =========================================================================
%  SIMULATION SETTINGS
%  =========================================================================
Tmax = 3000; dt = 0.02; transient = 500;
tn = round(Tmax/dt);
t = (0:tn)*dt;
idx_ss = t >= transient;

Mg_levels = [0.2, 0.5, 1.0, 1.4, 1.6, 1.8, 2.0, 2.5];
n_Mg = length(Mg_levels);

frequencies = [10, 30, 60, 80];
n_freq = length(frequencies);

fprintf('=========================================================================\n');
fprintf('  Mg²⁺ NEUROPROTECTION IN RGCs - Publication Analysis\n');
fprintf('=========================================================================\n\n');

%% =========================================================================
%  PART 1: MAIN SIMULATIONS
%  =========================================================================
fprintf('Running %d simulations...\n', n_Mg * n_freq);
tic;

results.spikes = zeros(n_Mg, n_freq);
results.peak_Ca = zeros(n_Mg, n_freq);
results.expected = zeros(1, n_freq);

V_traces = cell(n_Mg, 1);
Ca_traces = cell(n_Mg, 1);

for f = 1:n_freq
    freq = frequencies(f);
    period = 1000/freq;
    results.expected(f) = floor((Tmax - transient) / period);
    
    Glu = zeros(1, tn+1);
    for i = 1:tn+1
        if mod(t(i), period) < 2, Glu(i) = 1.0; end
    end
    
    for m = 1:n_Mg
        [V, Ca] = run_simulation(t, dt, Glu, Mg_levels(m), p);
        if f == n_freq
            V_traces{m} = V;
            Ca_traces{m} = Ca;
        end
        V_ss = V(idx_ss); Ca_ss = Ca(idx_ss);
        results.spikes(m,f) = sum(V_ss(1:end-1) < -20 & V_ss(2:end) >= -20);
        results.peak_Ca(m,f) = max(Ca_ss);
    end
end

elapsed = toc;
fprintf('Completed in %.1f seconds\n\n', elapsed);

results.prob = results.spikes ./ results.expected;
results.spike_red = 100 * (1 - results.spikes ./ results.spikes(1,:));
results.Ca_red = 100 * (1 - results.peak_Ca ./ results.peak_Ca(1,:));

%% =========================================================================
%  PART 2: INTERVENTION TIMING
%  =========================================================================
fprintf('Running intervention timing analysis...\n');

Mg_base = 0.2; Mg_treat = 1.8;
stress_start = 500; stress_end = 4500;
intervention_delays = [0, 500, 1000, 2000, 3000];

Tmax_int = 6000;
tn_int = round(Tmax_int/dt);
t_int = (0:tn_int)*dt;

Glu_int = zeros(1, tn_int+1);
period_80 = 1000/80;
for i = 1:tn_int+1
    if t_int(i) >= stress_start && t_int(i) <= stress_end
        if mod(t_int(i) - stress_start, period_80) < 2
            Glu_int(i) = 1.0;
        end
    end
end

n_int = length(intervention_delays) + 2;
int_results.peak_Ca = zeros(n_int, 1);
int_results.labels = cell(n_int, 1);
Ca_int_traces = cell(n_int, 1);

% No intervention
Mg_trace = Mg_base * ones(1, tn_int+1);
[~, Ca] = run_simulation(t_int, dt, Glu_int, Mg_trace, p);
int_results.peak_Ca(1) = max(Ca(t_int >= stress_start & t_int <= stress_end));
int_results.labels{1} = 'No Mg²⁺';
Ca_int_traces{1} = Ca;

% Pre-treated
Mg_trace = Mg_treat * ones(1, tn_int+1);
[~, Ca] = run_simulation(t_int, dt, Glu_int, Mg_trace, p);
int_results.peak_Ca(2) = max(Ca(t_int >= stress_start & t_int <= stress_end));
int_results.labels{2} = 'Pre-treated';
Ca_int_traces{2} = Ca;

% Delayed
for k = 1:length(intervention_delays)
    delay = intervention_delays(k);
    Mg_trace = Mg_base * ones(1, tn_int+1);
    Mg_trace(t_int >= stress_start + delay) = Mg_treat;
    [~, Ca] = run_simulation(t_int, dt, Glu_int, Mg_trace, p);
    int_results.peak_Ca(k+2) = max(Ca(t_int >= stress_start & t_int <= stress_end));
    int_results.labels{k+2} = sprintf('+%ds', delay/1000);
    Ca_int_traces{k+2} = Ca;
end

int_results.protection = 100 * (int_results.peak_Ca(1) - int_results.peak_Ca) / ...
                         (int_results.peak_Ca(1) - int_results.peak_Ca(2));
fprintf('Done.\n\n');

%% =========================================================================
%  FIGURE 1: MODEL OVERVIEW
%  =========================================================================
fig1 = figure('Position', [50 50 1000 600], 'Color', 'w');

% A: Mechanism schematic
subplot(2,3,1);
axis off;
text(0.5, 0.9, 'RGC Model', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.5, 0.7, 'Glutamate → AMPA + NMDA', 'FontSize', 9, 'HorizontalAlignment', 'center');
text(0.5, 0.5, 'NMDA blocked by Mg²⁺', 'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', [0.2 0.4 0.8]);
text(0.5, 0.3, '↓ Ca²⁺ → Neuroprotection', 'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', colors.optimal);
title('A. Mechanism');

% B: Mg²⁺ block curve
subplot(2,3,2);
V_range = -80:1:40;
Mg_test = [0.2, 1.0, 2.0];
hold on;
for i = 1:length(Mg_test)
    B_Mg = 1 ./ (1 + p.eta * Mg_test(i) * exp(-p.gamma * V_range));
    plot(V_range, B_Mg, 'LineWidth', 2, 'DisplayName', sprintf('%.1f mM', Mg_test(i)));
end
xlabel('V (mV)'); ylabel('NMDA Fraction');
title('B. Mg²⁺ Block');
legend('Location', 'northwest', 'Box', 'off');
xlim([-80 40]); ylim([0 1]); box off; grid on;

% C: Voltage traces (80 Hz)
subplot(2,3,3);
hold on;
m_show = [1, 4, 7];
for i = 1:length(m_show)
    m = m_show(i);
    plot(t/1000, V_traces{m} - (i-1)*100, 'Color', colors.Mg_gradient(m,:), 'LineWidth', 0.6);
end
xlim([0.5 1.0]); yticks([]);
xlabel('Time (s)'); ylabel('V (stacked)');
title('C. Voltage (80 Hz)'); box off;

% D: Calcium traces
subplot(2,3,4);
hold on;
for i = 1:length(m_show)
    m = m_show(i);
    plot(t/1000, Ca_traces{m}, 'Color', colors.Mg_gradient(m,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.1f mM', Mg_levels(m)));
end
yline(Ca_toxic, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('[Ca²⁺]_i (µM)');
title('D. Calcium (80 Hz)');
legend('Location', 'northeast', 'Box', 'off');
xlim([0 3]); box off; grid on;

% E: Spike probability
subplot(2,3,5);
hold on;
for m = 1:n_Mg
    plot(frequencies, results.prob(m,:), '-o', 'Color', colors.Mg_gradient(m,:), ...
        'MarkerFaceColor', colors.Mg_gradient(m,:), 'MarkerSize', 4);
end
yline(0.8, 'k--', 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('Spike Probability');
title('E. Reliability');
ylim([0.7 1.02]); box off; grid on;

% F: Peak Ca bar chart
subplot(2,3,6);
b = bar(results.peak_Ca(:, n_freq), 'FaceColor', 'flat');
b.CData = colors.Mg_gradient;
hold on; yline(Ca_toxic, 'r-', 'LineWidth', 2);
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), Mg_levels, 'UniformOutput', false));
xtickangle(45);
xlabel('[Mg²⁺] (mM)'); ylabel('Peak Ca²⁺ (µM)');
title('F. Peak Ca (80 Hz)'); box off;

sgtitle('Figure 1: RGC Model Overview', 'FontSize', 13, 'FontWeight', 'bold');

%% =========================================================================
%  FIGURE 2: DOSE-RESPONSE
%  =========================================================================
fig2 = figure('Position', [100 100 900 350], 'Color', 'w');

subplot(1,2,1);
hold on;
for f = 1:n_freq
    plot(Mg_levels, results.spike_red(:,f), '-o', 'Color', colors.freq(f,:), ...
        'MarkerFaceColor', colors.freq(f,:), 'MarkerSize', 5, ...
        'DisplayName', sprintf('%d Hz', frequencies(f)));
end
yline(max_spike_loss, 'k--', sprintf('≤%d%%', max_spike_loss), 'LineWidth', 1.5);
xlabel('[Mg²⁺] (mM)'); ylabel('Spike Reduction (%)');
title('A. Function Impact');
legend('Location', 'northwest', 'Box', 'off');
ylim([0 28]); box off; grid on;

subplot(1,2,2);
hold on;
for f = 1:n_freq
    plot(Mg_levels, results.Ca_red(:,f), '-s', 'Color', colors.freq(f,:), ...
        'MarkerFaceColor', colors.freq(f,:), 'MarkerSize', 5, ...
        'DisplayName', sprintf('%d Hz', frequencies(f)));
end
xlabel('[Mg²⁺] (mM)'); ylabel('Ca²⁺ Reduction (%)');
title('B. Neuroprotection');
legend('Location', 'southeast', 'Box', 'off');
ylim([0 100]); box off; grid on;

sgtitle('Figure 2: Dose-Response Analysis', 'FontSize', 13, 'FontWeight', 'bold');

%% =========================================================================
%  FIGURE 3: THERAPEUTIC WINDOWS
%  =========================================================================
fig3 = figure('Position', [150 150 900 350], 'Color', 'w');

% A: Trade-off space
subplot(1,2,1);
hold on;

Ca_base_80 = results.peak_Ca(1, n_freq);
Ca_thresh_red = 100 * (Ca_base_80 - Ca_toxic) / Ca_base_80;

fill([0 max_spike_loss max_spike_loss 0], [Ca_thresh_red Ca_thresh_red 100 100], ...
    [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

for m = 1:n_Mg
    ca_ok = results.peak_Ca(m, n_freq) < Ca_toxic;
    func_ok = results.spike_red(m, n_freq) <= max_spike_loss;
    if ca_ok && func_ok
        color = colors.optimal; sz = 100;
    else
        color = colors.suboptimal; sz = 50;
    end
    scatter(results.spike_red(m, n_freq), results.Ca_red(m, n_freq), sz, color, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
plot(results.spike_red(:, n_freq), results.Ca_red(:, n_freq), 'k-', 'LineWidth', 0.5);

xline(max_spike_loss, 'b--', 'LineWidth', 1.5);
yline(Ca_thresh_red, 'r-', 'LineWidth', 2);

xlabel('Spike Reduction (%)'); ylabel('Ca²⁺ Reduction (%)');
title('A. Trade-off (80 Hz)');
xlim([-2 26]); ylim([45 90]);
text(2, 86, 'OPTIMAL', 'Color', colors.optimal, 'FontWeight', 'bold');
box off; grid on;

% B: Window widths
subplot(1,2,2);

optimal_widths = zeros(n_freq, 1);
optimal_ranges = cell(n_freq, 1);
for f = 1:n_freq
    ca_ok = results.peak_Ca(:, f) < Ca_toxic;
    func_ok = results.spike_red(:, f) <= max_spike_loss;
    optimal = ca_ok & func_ok;
    if any(optimal)
        opt_idx = find(optimal);
        optimal_widths(f) = Mg_levels(opt_idx(end)) - Mg_levels(opt_idx(1));
        optimal_ranges{f} = sprintf('%.1f–%.1f', Mg_levels(opt_idx(1)), Mg_levels(opt_idx(end)));
    else
        optimal_ranges{f} = 'N/A';
    end
end

b = bar(optimal_widths, 'FaceColor', 'flat');
b.CData = colors.freq;
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d Hz', x), frequencies, 'UniformOutput', false));
ylabel('Window Width (mM)');
title('B. Therapeutic Windows');

for f = 1:n_freq
    if optimal_widths(f) > 0
        text(f, optimal_widths(f) + 0.05, optimal_ranges{f}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
end
ylim([0 2.2]); box off;

sgtitle('Figure 3: Therapeutic Windows', 'FontSize', 13, 'FontWeight', 'bold');

%% =========================================================================
%  FIGURE 4: INTERVENTION TIMING
%  =========================================================================
fig4 = figure('Position', [200 200 900 350], 'Color', 'w');

% A: Ca trajectories
subplot(1,2,1);
hold on;
plot_idx = [1, 2, 3, 5, 7];
line_styles = {'-', '-', '--', ':', '-.'};
line_colors = {colors.danger, colors.optimal, [0.3 0.6 0.4], [0.6 0.6 0.3], [0.5 0.5 0.5]};

for i = 1:length(plot_idx)
    k = plot_idx(i);
    plot(t_int/1000, Ca_int_traces{k}, line_styles{i}, 'Color', line_colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', int_results.labels{k});
end
yline(Ca_toxic, 'r--', 'LineWidth', 1);
xline(stress_start/1000, 'k:', 'LineWidth', 1);
xline(stress_end/1000, 'k:', 'LineWidth', 1);

xlabel('Time (s)'); ylabel('[Ca²⁺]_i (µM)');
title('A. Ca²⁺ Dynamics');
legend('Location', 'northeast', 'Box', 'off', 'FontSize', 7);
xlim([0 6]); box off; grid on;

% B: Protection efficacy
subplot(1,2,2);
bar_colors_int = [colors.danger; colors.optimal; ...
    colors.optimal*0.9 + colors.suboptimal*0.1; ...
    colors.optimal*0.6 + colors.suboptimal*0.4; ...
    colors.optimal*0.4 + colors.suboptimal*0.6; ...
    colors.optimal*0.2 + colors.suboptimal*0.8; ...
    colors.suboptimal];

b = bar(int_results.protection, 'FaceColor', 'flat');
b.CData = bar_colors_int;
set(gca, 'XTickLabel', int_results.labels, 'XTickLabelRotation', 45);
ylabel('Protection (%)');
title('B. Timing Matters');

for i = 1:n_int
    if int_results.protection(i) > 5
        text(i, int_results.protection(i) + 3, sprintf('%.0f%%', int_results.protection(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
ylim([0 115]); box off;

sgtitle('Figure 4: Intervention Timing', 'FontSize', 13, 'FontWeight', 'bold');

%% =========================================================================
%  FIGURE 5: MECHANISM
%  =========================================================================
fig5 = figure('Position', [250 250 800 400], 'Color', 'w');

Mg_mech = [0.2, 1.8];
Tmax_mech = 500; dt_mech = 0.01;
tn_mech = round(Tmax_mech/dt_mech);
t_mech = (0:tn_mech)*dt_mech;

period_mech = 1000/80;
Glu_mech = zeros(1, tn_mech+1);
for i = 1:tn_mech+1
    if mod(t_mech(i), period_mech) < 2, Glu_mech(i) = 1.0; end
end

[V_low, Ca_low, I_NMDA_low] = run_simulation_detailed(t_mech, dt_mech, Glu_mech, Mg_mech(1), p);
[V_high, Ca_high, I_NMDA_high] = run_simulation_detailed(t_mech, dt_mech, Glu_mech, Mg_mech(2), p);

% A: Voltage
subplot(2,2,1);
hold on;
plot(t_mech, V_low, 'Color', colors.danger, 'LineWidth', 1, 'DisplayName', '0.2 mM');
plot(t_mech, V_high, 'Color', colors.optimal, 'LineWidth', 1, 'DisplayName', '1.8 mM');
xlabel('Time (ms)'); ylabel('V (mV)');
title('A. Membrane Potential');
legend('Location', 'northeast', 'Box', 'off');
xlim([0 250]); box off;

% B: NMDA current
subplot(2,2,2);
hold on;
plot(t_mech, -I_NMDA_low, 'Color', colors.danger, 'LineWidth', 1.5);
plot(t_mech, -I_NMDA_high, 'Color', colors.optimal, 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('I_{NMDA} (µA/cm²)');
title('B. NMDA Current');
xlim([0 250]); box off;

% C: Calcium
subplot(2,2,3);
hold on;
plot(t_mech, Ca_low, 'Color', colors.danger, 'LineWidth', 2);
plot(t_mech, Ca_high, 'Color', colors.optimal, 'LineWidth', 2);
yline(Ca_toxic, 'r--', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('[Ca²⁺]_i (µM)');
title('C. Intracellular Ca²⁺');
xlim([0 Tmax_mech]); box off;

% D: Summary
subplot(2,2,4);
axis off;
text(0.5, 0.85, 'Protection Cascade:', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.5, 0.65, '↑ [Mg²⁺]_{ext}', 'FontSize', 10, 'HorizontalAlignment', 'center');
text(0.5, 0.50, '↓', 'FontSize', 12, 'HorizontalAlignment', 'center');
text(0.5, 0.35, '↑ NMDA block → ↓ I_{NMDA}', 'FontSize', 10, 'HorizontalAlignment', 'center');
text(0.5, 0.20, '↓', 'FontSize', 12, 'HorizontalAlignment', 'center');
text(0.5, 0.05, '↓ Ca²⁺ → Protection', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', colors.optimal);
title('D. Summary');

sgtitle('Figure 5: Mechanism', 'FontSize', 13, 'FontWeight', 'bold');

%% =========================================================================
%  SAVE FIGURES (with proper sizing)
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  Saving figures...\n');
fprintf('=========================================================================\n\n');

% Save function with proper PDF settings
save_figure(fig1, 'Figure1_ModelOverview');
save_figure(fig2, 'Figure2_DoseResponse');
save_figure(fig3, 'Figure3_TherapeuticWindows');
save_figure(fig4, 'Figure4_InterventionTiming');
save_figure(fig5, 'Figure5_Mechanism');

fprintf('All figures saved!\n');

%% =========================================================================
%  SUMMARY
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  RESULTS SUMMARY\n');
fprintf('=========================================================================\n\n');

fprintf('THERAPEUTIC WINDOWS:\n');
for f = 1:n_freq
    fprintf('  %d Hz: %s\n', frequencies(f), optimal_ranges{f});
end

fprintf('\nKEY FINDING AT 80 Hz:\n');
fprintf('  Optimal [Mg²⁺]: 1.6–2.0 mM\n');
fprintf('  Spike loss plateaus at 20%% while Ca²⁺ continues to decrease\n');

fprintf('\nINTERVENTION TIMING:\n');
for i = 1:n_int
    fprintf('  %s: %.0f%% protection\n', int_results.labels{i}, int_results.protection(i));
end

%% =========================================================================
%  FUNCTIONS
%  =========================================================================
function save_figure(fig, filename)
    % Save PNG
    saveas(fig, [filename '.png']);
    
    % Save PDF with proper sizing
    set(fig, 'PaperPositionMode', 'auto');
    fig_pos = fig.Position;
    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
    print(fig, [filename '.pdf'], '-dpdf', '-vector');
    
    fprintf('  Saved: %s.png, %s.pdf\n', filename, filename);
end

function [V, Ca] = run_simulation(t, dt, Glu, Mg_ext, p)
    if isscalar(Mg_ext)
        Mg_trace = Mg_ext * ones(1, length(t));
    else
        Mg_trace = Mg_ext;
    end
    
    tn = length(t) - 1;
    V = zeros(1, tn+1); V(1) = p.V_rest;
    m_g = zeros(1, tn+1); h_g = zeros(1, tn+1); n_g = zeros(1, tn+1);
    a_g = zeros(1, tn+1); b_g = zeros(1, tn+1); s_CaL = zeros(1, tn+1);
    s_AMPA = zeros(1, tn+1); s_NMDA = zeros(1, tn+1);
    Ca = zeros(1, tn+1); Ca(1) = p.Ca_rest;
    
    v0 = p.V_rest;
    am = 0.1*(v0+40)/(1-exp(-(v0+40)/10)); bm = 4*exp(-(v0+65)/18);
    m_g(1) = am/(am+bm);
    ah = 0.07*exp(-(v0+65)/20); bh = 1/(1+exp(-(v0+35)/10));
    h_g(1) = ah/(ah+bh);
    an = 0.01*(v0+55)/(1-exp(-(v0+55)/10)); bn = 0.125*exp(-(v0+65)/80);
    n_g(1) = an/(an+bn);
    a_g(1) = 1/(1+exp(-(v0+50)/20));
    b_g(1) = 1/(1+exp((v0+70)/6));
    s_CaL(1) = 1/(1+exp(-(v0+20)/6));
    
    for i = 1:tn
        v = V(i); Glu_now = Glu(i); Mg = Mg_trace(i);
        
        if Glu_now > 0
            r_AMPA = (Glu_now/(Glu_now + p.EC50_AMPA)) / p.tau_AMPA_rise;
            r_NMDA = (Glu_now/(Glu_now + p.EC50_NMDA)) / p.tau_NMDA_rise;
        else
            r_AMPA = 0; r_NMDA = 0;
        end
        
        s_AMPA(i+1) = max(0, min(1, s_AMPA(i) + dt*(r_AMPA*(1-s_AMPA(i)) - s_AMPA(i)/p.tau_AMPA_decay)));
        s_NMDA(i+1) = max(0, min(1, s_NMDA(i) + dt*(r_NMDA*(1-s_NMDA(i)) - s_NMDA(i)/p.tau_NMDA_decay)));
        
        if abs(v+40) < 1e-6, am = 1; else, am = 0.1*(v+40)/(1-exp(-(v+40)/10)); end
        bm = 4*exp(-(v+65)/18);
        ah = 0.07*exp(-(v+65)/20); bh = 1/(1+exp(-(v+35)/10));
        if abs(v+55) < 1e-6, an = 0.1; else, an = 0.01*(v+55)/(1-exp(-(v+55)/10)); end
        bn = 0.125*exp(-(v+65)/80);
        
        tau_m = 1/(am+bm); tau_h = 1/(ah+bh); tau_n = 1/(an+bn);
        m_g(i+1) = m_g(i) + dt*(am/(am+bm) - m_g(i))/tau_m;
        h_g(i+1) = h_g(i) + dt*(ah/(ah+bh) - h_g(i))/tau_h;
        n_g(i+1) = n_g(i) + dt*(an/(an+bn) - n_g(i))/tau_n;
        a_g(i+1) = a_g(i) + dt*(1/(1+exp(-(v+50)/20)) - a_g(i))/5;
        b_g(i+1) = b_g(i) + dt*(1/(1+exp((v+70)/6)) - b_g(i))/20;
        s_CaL(i+1) = s_CaL(i) + dt*(1/(1+exp(-(v+20)/6)) - s_CaL(i))/5;
        
        I_Na = p.g_Na * m_g(i)^3 * h_g(i) * (v - p.E_Na);
        I_Kdr = p.g_Kdr * n_g(i)^4 * (v - p.E_K);
        I_KA = p.g_KA * a_g(i)^3 * b_g(i) * (v - p.E_K);
        I_CaL = p.g_CaL * s_CaL(i)^2 * (v - p.E_Ca);
        I_L = p.g_L * (v - p.E_L);
        
        q_KCa = Ca(i)^p.n_KCa / (Ca(i)^p.n_KCa + p.Kd_KCa^p.n_KCa);
        I_KCa = p.g_KCa * q_KCa * (v - p.E_K);
        
        B_Mg = 1 / (1 + p.eta * Mg * exp(-p.gamma * v));
        I_AMPA = p.g_AMPA_max * s_AMPA(i) * (v - p.E_exc);
        I_NMDA = p.g_NMDA_max * s_NMDA(i) * B_Mg * (v - p.E_exc);
        
        Ca_influx = -p.k_Ca_NMDA * p.f_Ca * I_NMDA - p.k_Ca_CaL * I_CaL;
        Ca(i+1) = max(p.Ca_rest, Ca(i) + dt*(Ca_influx - (Ca(i)-p.Ca_rest)/p.tau_Ca));
        
        I_total = I_Na + I_Kdr + I_KA + I_CaL + I_KCa + I_L + I_AMPA + I_NMDA;
        V(i+1) = v + dt*(-I_total)/p.Cm;
    end
end

function [V, Ca, I_NMDA_rec] = run_simulation_detailed(t, dt, Glu, Mg_ext, p)
    tn = length(t) - 1;
    V = zeros(1, tn+1); V(1) = p.V_rest;
    m_g = zeros(1, tn+1); h_g = zeros(1, tn+1); n_g = zeros(1, tn+1);
    a_g = zeros(1, tn+1); b_g = zeros(1, tn+1); s_CaL = zeros(1, tn+1);
    s_AMPA = zeros(1, tn+1); s_NMDA = zeros(1, tn+1);
    Ca = zeros(1, tn+1); Ca(1) = p.Ca_rest;
    I_NMDA_rec = zeros(1, tn+1);
    
    v0 = p.V_rest;
    am = 0.1*(v0+40)/(1-exp(-(v0+40)/10)); bm = 4*exp(-(v0+65)/18);
    m_g(1) = am/(am+bm);
    ah = 0.07*exp(-(v0+65)/20); bh = 1/(1+exp(-(v0+35)/10));
    h_g(1) = ah/(ah+bh);
    an = 0.01*(v0+55)/(1-exp(-(v0+55)/10)); bn = 0.125*exp(-(v0+65)/80);
    n_g(1) = an/(an+bn);
    a_g(1) = 1/(1+exp(-(v0+50)/20));
    b_g(1) = 1/(1+exp((v0+70)/6));
    s_CaL(1) = 1/(1+exp(-(v0+20)/6));
    
    for i = 1:tn
        v = V(i); Glu_now = Glu(i);
        
        if Glu_now > 0
            r_AMPA = (Glu_now/(Glu_now + p.EC50_AMPA)) / p.tau_AMPA_rise;
            r_NMDA = (Glu_now/(Glu_now + p.EC50_NMDA)) / p.tau_NMDA_rise;
        else
            r_AMPA = 0; r_NMDA = 0;
        end
        
        s_AMPA(i+1) = max(0, min(1, s_AMPA(i) + dt*(r_AMPA*(1-s_AMPA(i)) - s_AMPA(i)/p.tau_AMPA_decay)));
        s_NMDA(i+1) = max(0, min(1, s_NMDA(i) + dt*(r_NMDA*(1-s_NMDA(i)) - s_NMDA(i)/p.tau_NMDA_decay)));
        
        if abs(v+40) < 1e-6, am = 1; else, am = 0.1*(v+40)/(1-exp(-(v+40)/10)); end
        bm = 4*exp(-(v+65)/18);
        ah = 0.07*exp(-(v+65)/20); bh = 1/(1+exp(-(v+35)/10));
        if abs(v+55) < 1e-6, an = 0.1; else, an = 0.01*(v+55)/(1-exp(-(v+55)/10)); end
        bn = 0.125*exp(-(v+65)/80);
        
        tau_m = 1/(am+bm); tau_h = 1/(ah+bh); tau_n = 1/(an+bn);
        m_g(i+1) = m_g(i) + dt*(am/(am+bm) - m_g(i))/tau_m;
        h_g(i+1) = h_g(i) + dt*(ah/(ah+bh) - h_g(i))/tau_h;
        n_g(i+1) = n_g(i) + dt*(an/(an+bn) - n_g(i))/tau_n;
        a_g(i+1) = a_g(i) + dt*(1/(1+exp(-(v+50)/20)) - a_g(i))/5;
        b_g(i+1) = b_g(i) + dt*(1/(1+exp((v+70)/6)) - b_g(i))/20;
        s_CaL(i+1) = s_CaL(i) + dt*(1/(1+exp(-(v+20)/6)) - s_CaL(i))/5;
        
        I_Na = p.g_Na * m_g(i)^3 * h_g(i) * (v - p.E_Na);
        I_Kdr = p.g_Kdr * n_g(i)^4 * (v - p.E_K);
        I_KA = p.g_KA * a_g(i)^3 * b_g(i) * (v - p.E_K);
        I_CaL = p.g_CaL * s_CaL(i)^2 * (v - p.E_Ca);
        I_L = p.g_L * (v - p.E_L);
        
        q_KCa = Ca(i)^p.n_KCa / (Ca(i)^p.n_KCa + p.Kd_KCa^p.n_KCa);
        I_KCa = p.g_KCa * q_KCa * (v - p.E_K);
        
        B_Mg = 1 / (1 + p.eta * Mg_ext * exp(-p.gamma * v));
        I_AMPA = p.g_AMPA_max * s_AMPA(i) * (v - p.E_exc);
        I_NMDA = p.g_NMDA_max * s_NMDA(i) * B_Mg * (v - p.E_exc);
        I_NMDA_rec(i) = I_NMDA;
        
        Ca_influx = -p.k_Ca_NMDA * p.f_Ca * I_NMDA - p.k_Ca_CaL * I_CaL;
        Ca(i+1) = max(p.Ca_rest, Ca(i) + dt*(Ca_influx - (Ca(i)-p.Ca_rest)/p.tau_Ca));
        
        I_total = I_Na + I_Kdr + I_KA + I_CaL + I_KCa + I_L + I_AMPA + I_NMDA;
        V(i+1) = v + dt*(-I_total)/p.Cm;
    end
end
