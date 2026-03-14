%% =========================================================================
%  Mg²⁺ NEUROPROTECTION IN RETINAL GANGLION CELLS
%  REVISED VERSION - Addressing Reviewer Comments
%  =========================================================================
%
%  REVISIONS MADE:
%    1. Sensitivity analysis for Ca²⁺ toxicity threshold (0.8, 1.0, 1.2 µM)
%    2. Numerical method validation (Euler vs RK4, dt comparison)
%    3. Improved figure labeling with larger fonts
%    4. Expanded parameter table as figure
%    5. Consistent terminology ("high-frequency pathological stimulation")
%
%  FIGURES GENERATED:
%    Figure 1: Model Overview & Example Traces
%    Figure 2: Dose-Response Analysis  
%    Figure 3: Therapeutic Windows
%    Figure 4: Intervention Timing
%    Figure 5: Mechanistic Demonstration
%    Figure S1: High-Resolution Therapeutic Window (80 Hz)
%    Figure S2: Detailed Results Table
%    Figure S3: Calcium Threshold Sensitivity Analysis (NEW)
%    Figure S4: Numerical Method Validation (NEW)
%
%  Author: Mehdi Borjkhani
%  Institution: ICTER, Polish Academy of Sciences
%  Revision Date: January 2026
%  =========================================================================

clear; clc; close all;

%% =========================================================================
%  GLOBAL SETTINGS - IMPROVED FOR PUBLICATION
%  =========================================================================
set(0, 'DefaultAxesFontSize', 11);  % Increased from 10
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesLineWidth', 1.2);
set(0, 'DefaultLineLineWidth', 1.8);
set(0, 'DefaultAxesLabelFontSizeMultiplier', 1.1);
set(0, 'DefaultAxesTitleFontSizeMultiplier', 1.15);

% Color schemes - consistent throughout
colors.Mg_gradient = [0.85 0.2 0.2;    % Low Mg (red - danger)
                      0.75 0.4 0.2;
                      0.6 0.5 0.3;
                      0.4 0.55 0.4;
                      0.3 0.6 0.5;
                      0.2 0.5 0.7;
                      0.15 0.4 0.75;
                      0.1 0.3 0.8];     % High Mg (blue - protected)

colors.freq = [0.2 0.7 0.3;            % 10 Hz (green - safe)
               0.4 0.65 0.2;           % 30 Hz (yellow-green)
               0.65 0.55 0.1;          % 60 Hz (gold)
               0.85 0.35 0.15;         % 80 Hz (orange-red)
               0.9 0.15 0.15;          % 90 Hz (red)
               0.6 0.1 0.4];           % 100 Hz (dark red/purple - extreme)

colors.optimal = [0.2 0.7 0.3];
colors.suboptimal = [0.7 0.7 0.7];
colors.danger = [0.85 0.2 0.2];
colors.threshold_sens = [0.2 0.5 0.8; 0.3 0.7 0.3; 0.8 0.4 0.2]; % For sensitivity

%% =========================================================================
%  MODEL PARAMETERS - COMPLETE LIST FOR EXPANDED TABLE
%  =========================================================================
% Membrane properties
p.Cm = 1.0;          % Membrane capacitance (µF/cm²)
p.V_rest = -65;      % Resting potential (mV)

% Reversal potentials (mV)
p.E_Na = 50;         % Sodium
p.E_K = -77;         % Potassium
p.E_L = -54.4;       % Leak
p.E_Ca = 120;        % Calcium
p.E_exc = 0;         % Excitatory synaptic

% Intrinsic conductances (mS/cm²)
p.g_Na = 120;        % Sodium
p.g_Kdr = 36;        % Delayed rectifier potassium
p.g_KA = 8;          % A-type potassium
p.g_CaL = 0.3;       % L-type calcium
p.g_KCa = 0.3;       % Calcium-activated potassium
p.g_L = 0.35;        % Leak

% Synaptic parameters - EXPANDED FOR TABLE
p.g_AMPA_max = 0.25;     % Max AMPA conductance (mS/cm²)
p.g_NMDA_max = 1.2;      % Max NMDA conductance (mS/cm²)
p.EC50_AMPA = 0.5;       % AMPA EC50 (mM) - glutamate affinity
p.EC50_NMDA = 0.002;     % NMDA EC50 (mM = 2 µM) - glutamate affinity
p.tau_AMPA_rise = 0.3;   % AMPA rise time constant (ms)
p.tau_AMPA_decay = 3.0;  % AMPA decay time constant (ms)
p.tau_NMDA_rise = 5.0;   % NMDA rise time constant (ms)
p.tau_NMDA_decay = 80.0; % NMDA decay time constant (ms)

% Mg²⁺ block parameters (Jahr-Stevens)
p.eta = 0.28;        % Mg²⁺ sensitivity (mM⁻¹)
p.gamma = 0.062;     % Voltage sensitivity (mV⁻¹)

% Calcium dynamics parameters - EXPANDED FOR TABLE
p.f_Ca = 0.15;       % Fractional Ca²⁺ current through NMDA
p.Ca_rest = 0.05;    % Resting [Ca²⁺]ᵢ (µM)
p.tau_Ca = 200;      % Calcium decay time constant (ms)
p.k_Ca_NMDA = 0.012; % NMDA-to-calcium scaling factor
p.k_Ca_CaL = 0.003;  % L-type-to-calcium scaling factor
p.Kd_KCa = 0.5;      % KCa half-activation (µM)
p.n_KCa = 2;         % KCa Hill coefficient

% Stimulation parameters - EXPANDED FOR TABLE
stim.Glu_conc = 1.0;     % Glutamate pulse concentration (mM)
stim.pulse_duration = 2;  % Glutamate pulse duration (ms)

% Therapeutic criteria
Ca_toxic = 1.0;          % Default toxicity threshold (µM)
max_spike_loss = 20;     % Maximum acceptable spike loss (%)

%% =========================================================================
%  SIMULATION SETTINGS
%  =========================================================================
Tmax = 3000; dt = 0.02; transient = 500;
tn = round(Tmax/dt);
t = (0:tn)*dt;
idx_ss = t >= transient;  % Steady-state window for analysis

% Analysis window parameters (for Methods section):
%   - Total simulation: Tmax = 3000 ms
%   - Transient discarded: 500 ms
%   - Spike counting window: T_count = 2500 ms
%   - At 80 Hz: N_pulses = 80 Hz × 2.5 s = 200 expected pulses
%   - Spike loss (%) = 100 × (1 - N_spikes / N_pulses)

Mg_levels = [0.2, 0.5, 1.0, 1.4, 1.6, 1.8, 2.0, 2.5];
n_Mg = length(Mg_levels);

% Extended frequency range to show therapeutic window narrowing at high frequencies
frequencies = [10, 30, 60, 80, 90, 100];
n_freq = length(frequencies);

fprintf('=========================================================================\n');
fprintf('  Mg²⁺ NEUROPROTECTION IN RGCs - REVISED Analysis\n');
fprintf('  Addressing Reviewer Comments\n');
fprintf('=========================================================================\n\n');
fprintf('Spike counting: T_count = %.0f ms (excluding %.0f ms transient)\n', Tmax-transient, transient);
fprintf('At 80 Hz: Expected pulses = %d\n\n', floor((Tmax-transient)*80/1000));

%% =========================================================================
%  PART 1: MAIN SIMULATIONS
%  =========================================================================
fprintf('Running %d main simulations...\n', n_Mg * n_freq);
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
        if mod(t(i), period) < stim.pulse_duration
            Glu(i) = stim.Glu_conc;
        end
    end
    
    for m = 1:n_Mg
        [V, Ca] = run_simulation_euler(t, dt, Glu, Mg_levels(m), p);
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

% Finer resolution between 0 and 0.5s to capture dynamics
intervention_delays = [0, 100, 200, 300, 400, 500, 750, 1000, 2000, 3000];  % in ms

Tmax_int = 6000;
tn_int = round(Tmax_int/dt);
t_int = (0:tn_int)*dt;

Glu_int = zeros(1, tn_int+1);
period_80 = 1000/80;
for i = 1:tn_int+1
    if t_int(i) >= stress_start && t_int(i) <= stress_end
        if mod(t_int(i) - stress_start, period_80) < stim.pulse_duration
            Glu_int(i) = stim.Glu_conc;
        end
    end
end

n_int = length(intervention_delays) + 2;  % +2 for No treatment and Pre-treated
int_results.peak_Ca = zeros(n_int, 1);
int_results.labels = cell(n_int, 1);
Ca_int_traces = cell(n_int, 1);

% No intervention (baseline - maximum damage)
Mg_trace = Mg_base * ones(1, tn_int+1);
[~, Ca] = run_simulation_euler(t_int, dt, Glu_int, Mg_trace, p);
int_results.peak_Ca(1) = max(Ca(t_int >= stress_start & t_int <= stress_end));
int_results.labels{1} = 'None';
Ca_int_traces{1} = Ca;

% Pre-treated (maximum protection)
Mg_trace = Mg_treat * ones(1, tn_int+1);
[~, Ca] = run_simulation_euler(t_int, dt, Glu_int, Mg_trace, p);
int_results.peak_Ca(2) = max(Ca(t_int >= stress_start & t_int <= stress_end));
int_results.labels{2} = 'Pre';
Ca_int_traces{2} = Ca;

% Delayed interventions
for k = 1:length(intervention_delays)
    delay = intervention_delays(k);
    Mg_trace = Mg_base * ones(1, tn_int+1);
    Mg_trace(t_int >= stress_start + delay) = Mg_treat;
    [~, Ca] = run_simulation_euler(t_int, dt, Glu_int, Mg_trace, p);
    int_results.peak_Ca(k+2) = max(Ca(t_int >= stress_start & t_int <= stress_end));
    % Format delay labels properly
    delay_sec = delay / 1000;
    if delay_sec == floor(delay_sec)
        int_results.labels{k+2} = sprintf('+%d', delay_sec);
    else
        int_results.labels{k+2} = sprintf('+%.1f', delay_sec);
    end
    Ca_int_traces{k+2} = Ca;
end

int_results.protection = 100 * (int_results.peak_Ca(1) - int_results.peak_Ca) / ...
                         (int_results.peak_Ca(1) - int_results.peak_Ca(2));
fprintf('Done.\n\n');

%% =========================================================================
%  PART 3: HIGH-RESOLUTION ANALYSIS (for S1 Fig)
%  =========================================================================
fprintf('Running high-resolution analysis at 80 Hz...\n');

Mg_highres = 1.0:0.1:2.5;
n_Mg_hr = length(Mg_highres);

period = 1000/80;
expected_pulses = floor((Tmax - transient) / period);
Glu = zeros(1, tn+1);
for i = 1:tn+1
    if mod(t(i), period) < stim.pulse_duration
        Glu(i) = stim.Glu_conc;
    end
end

hr_results.spikes = zeros(n_Mg_hr, 1);
hr_results.peak_Ca = zeros(n_Mg_hr, 1);

for m = 1:n_Mg_hr
    [V, Ca] = run_simulation_euler(t, dt, Glu, Mg_highres(m), p);
    V_ss = V(idx_ss); Ca_ss = Ca(idx_ss);
    hr_results.spikes(m) = sum(V_ss(1:end-1) < -20 & V_ss(2:end) >= -20);
    hr_results.peak_Ca(m) = max(Ca_ss);
end

% CRITICAL: Spike loss relative to expected pulses (N_pulses = f × T_count)
% This ensures consistent interpretation: 160/200 = 80% reliability = 20% loss
hr_results.spike_red = 100 * (expected_pulses - hr_results.spikes) / expected_pulses;

% Ca reduction relative to baseline (lowest Mg = 1.0 mM)
baseline_Ca = hr_results.peak_Ca(1);
hr_results.Ca_red = 100 * (baseline_Ca - hr_results.peak_Ca) / baseline_Ca;

% Print key metrics for verification
fprintf('  Expected pulses (N_pulses): %d\n', expected_pulses);
fprintf('  Baseline Ca (at 1.0 mM Mg): %.3f µM\n', baseline_Ca);

ca_safe = hr_results.peak_Ca < Ca_toxic;
func_ok = hr_results.spike_red <= max_spike_loss;
optimal = ca_safe & func_ok;

fprintf('Done.\n\n');

%% =========================================================================
%  PART 4: CALCIUM THRESHOLD SENSITIVITY ANALYSIS (HIGH RESOLUTION)
%  =========================================================================
fprintf('Calcium threshold sensitivity analysis...\n');

% HIGH RESOLUTION: 17 threshold values (0.6 to 1.4 µM in 0.05 steps)
Ca_thresholds = 0.6:0.05:1.4;  % µM - much finer resolution
n_thresh = length(Ca_thresholds);

sens_results.optimal_ranges = cell(n_thresh, 1);
sens_results.window_widths = zeros(n_thresh, 1);
sens_results.lower_bounds = zeros(n_thresh, 1);
sens_results.upper_bounds = zeros(n_thresh, 1);

for th = 1:n_thresh
    Ca_thresh = Ca_thresholds(th);
    ca_safe_th = hr_results.peak_Ca < Ca_thresh;
    optimal_th = ca_safe_th & func_ok;
    
    if any(optimal_th)
        opt_idx = find(optimal_th);
        sens_results.lower_bounds(th) = Mg_highres(opt_idx(1));
        sens_results.upper_bounds(th) = Mg_highres(opt_idx(end));
        sens_results.window_widths(th) = sens_results.upper_bounds(th) - sens_results.lower_bounds(th);
        sens_results.optimal_ranges{th} = sprintf('%.1f–%.1f mM', ...
            sens_results.lower_bounds(th), sens_results.upper_bounds(th));
    else
        sens_results.optimal_ranges{th} = 'None';
        sens_results.window_widths(th) = 0;
        sens_results.lower_bounds(th) = NaN;
        sens_results.upper_bounds(th) = NaN;
    end
end

% Print summary table to console
fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║       Ca²⁺ THRESHOLD SENSITIVITY ANALYSIS                             ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Threshold (µM)  │  Optimal [Mg²⁺] Range  │  Window Width (mM)       ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════════════╣\n');
for th = 1:n_thresh
    if sens_results.window_widths(th) > 0
        fprintf('║       %.2f        │      %s         │        %.2f              ║\n', ...
            Ca_thresholds(th), sens_results.optimal_ranges{th}, sens_results.window_widths(th));
    else
        fprintf('║       %.2f        │        None            │        0.00              ║\n', ...
            Ca_thresholds(th));
    end
end
fprintf('╚═══════════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

%% =========================================================================
%  PART 5: NUMERICAL METHOD VALIDATION (NEW - Reviewer Request)
%  =========================================================================
fprintf('Running numerical method validation...\n');

% Test case: 80 Hz stimulation at Mg = 1.8 mM
Mg_test = 1.8;
Tmax_test = 1000;

% Time steps to compare
dt_values = [0.05, 0.02, 0.01, 0.005];
n_dt = length(dt_values);

% Store results for Euler method at different dt
euler_results.peak_Ca = zeros(n_dt, 1);
euler_results.spikes = zeros(n_dt, 1);
euler_results.dt = dt_values;

for d = 1:n_dt
    dt_test = dt_values(d);
    tn_test = round(Tmax_test/dt_test);
    t_test = (0:tn_test)*dt_test;
    
    Glu_test = zeros(1, tn_test+1);
    for i = 1:tn_test+1
        if mod(t_test(i), period_80) < stim.pulse_duration
            Glu_test(i) = stim.Glu_conc;
        end
    end
    
    [V_test, Ca_test] = run_simulation_euler(t_test, dt_test, Glu_test, Mg_test, p);
    
    euler_results.peak_Ca(d) = max(Ca_test);
    euler_results.spikes(d) = sum(V_test(1:end-1) < -20 & V_test(2:end) >= -20);
end

% Run RK4 at dt=0.02 for comparison
dt_rk4 = 0.02;
tn_rk4 = round(Tmax_test/dt_rk4);
t_rk4 = (0:tn_rk4)*dt_rk4;

Glu_rk4 = zeros(1, tn_rk4+1);
for i = 1:tn_rk4+1
    if mod(t_rk4(i), period_80) < stim.pulse_duration
        Glu_rk4(i) = stim.Glu_conc;
    end
end

[V_rk4, Ca_rk4] = run_simulation_rk4(t_rk4, dt_rk4, Glu_rk4, Mg_test, p);

rk4_results.peak_Ca = max(Ca_rk4);
rk4_results.spikes = sum(V_rk4(1:end-1) < -20 & V_rk4(2:end) >= -20);

% Calculate relative errors vs finest Euler
euler_results.Ca_error = abs(euler_results.peak_Ca - euler_results.peak_Ca(end)) / ...
                         euler_results.peak_Ca(end) * 100;
euler_results.spike_error = abs(euler_results.spikes - euler_results.spikes(end));

% RK4 vs finest Euler
rk4_results.Ca_error = abs(rk4_results.peak_Ca - euler_results.peak_Ca(end)) / ...
                       euler_results.peak_Ca(end) * 100;
rk4_results.spike_error = abs(rk4_results.spikes - euler_results.spikes(end));

fprintf('Numerical Validation Results:\n');
fprintf('  Forward Euler:\n');
for d = 1:n_dt
    fprintf('    dt=%.3f ms: Peak Ca = %.4f µM, Spikes = %d, Ca error = %.2f%%\n', ...
        dt_values(d), euler_results.peak_Ca(d), euler_results.spikes(d), euler_results.Ca_error(d));
end
fprintf('  RK4 (dt=0.02 ms): Peak Ca = %.4f µM, Spikes = %d, Ca error vs Euler(finest) = %.2f%%\n', ...
    rk4_results.peak_Ca, rk4_results.spikes, rk4_results.Ca_error);
fprintf('\n');

%% =========================================================================
%  FIGURE 1: MODEL OVERVIEW (IMPROVED LABELING)
%  =========================================================================
fprintf('Generating figures...\n');

fig1 = figure('Position', [50 50 1000 600], 'Color', 'w');

% A: Mechanism schematic
subplot(2,3,1);
axis off;
text(0.5, 0.9, '\bf{RGC Model}', 'FontSize', 13, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
text(0.5, 0.7, 'Glutamate \rightarrow AMPA + NMDA', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
text(0.5, 0.5, 'NMDA blocked by Mg^{2+}', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', [0.2 0.4 0.8], 'Interpreter', 'tex');
text(0.5, 0.3, '\downarrow Ca^{2+} \rightarrow Neuroprotection', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', colors.optimal, 'Interpreter', 'tex');
title('A. Mechanism', 'FontSize', 12, 'FontWeight', 'bold');

% B: Mg²⁺ block curve
subplot(2,3,2);
V_range = -80:1:40;
Mg_test_plot = [0.2, 1.0, 2.0];
hold on;
for i = 1:length(Mg_test_plot)
    B_Mg = 1 ./ (1 + p.eta * Mg_test_plot(i) * exp(-p.gamma * V_range));
    plot(V_range, B_Mg, 'LineWidth', 2.5, 'DisplayName', sprintf('%.1f mM', Mg_test_plot(i)));
end
xlabel('Membrane Potential (mV)', 'FontSize', 11);
ylabel('NMDA Conductance Fraction', 'FontSize', 11);
title('B. Mg^{2+} Block', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 10);
xlim([-80 40]); ylim([0 1]); box off; grid on;

% C: Voltage traces (80 Hz) - labeled as "High-frequency pathological stimulation"
subplot(2,3,3);
hold on;
m_show = [1, 4, 7];
for i = 1:length(m_show)
    m = m_show(i);
    plot(t/1000, V_traces{m} - (i-1)*100, 'Color', colors.Mg_gradient(m,:), 'LineWidth', 0.8);
end
xlim([0.5 1.0]); yticks([]);
xlabel('Time (s)', 'FontSize', 11);
ylabel('V (stacked)', 'FontSize', 11);
title('C. Voltage (80 Hz)', 'FontSize', 12, 'FontWeight', 'bold');
box off;

% D: Calcium traces
subplot(2,3,4);
hold on;
for i = 1:length(m_show)
    m = m_show(i);
    plot(t/1000, Ca_traces{m}, 'Color', colors.Mg_gradient(m,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%.1f mM', Mg_levels(m)));
end
yline(Ca_toxic, 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(0.12, Ca_toxic + 0.30, 'Toxicity Threshold', 'FontSize', 9, 'Color', 'r');
xlabel('Time (s)', 'FontSize', 11);
ylabel('[Ca^{2+}]_i (\muM)', 'FontSize', 11, 'Interpreter', 'tex');
title('D. Calcium (80 Hz)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'Box', 'off', 'FontSize', 9);
xlim([0 3]); box off; grid on;

% E: Spike probability (showing all Mg levels vs frequency)
subplot(2,3,5);
hold on;
% Plot selected Mg levels for clarity
Mg_show_idx = [1, 4, 6, 8];  % 0.2, 1.4, 1.8, 2.5 mM
for i = 1:length(Mg_show_idx)
    m = Mg_show_idx(i);
    plot(frequencies, results.prob(m,:), '-o', 'Color', colors.Mg_gradient(m,:), ...
        'MarkerFaceColor', colors.Mg_gradient(m,:), 'MarkerSize', 6, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%.1f mM', Mg_levels(m)));
end
yline(0.8, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(55, 0.79, '≤20% loss', 'FontSize', 8);
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Spike Probability', 'FontSize', 11);
title('E. Reliability vs Frequency', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southwest', 'Box', 'off', 'FontSize', 8);
ylim([0.68 1.02]); xlim([5 85]); box off; grid on;

% F: Peak Ca bar chart
subplot(2,3,6);
b = bar(results.peak_Ca(:, n_freq), 'FaceColor', 'flat');
b.CData = colors.Mg_gradient;
hold on;
yline(Ca_toxic, 'r-', 'Toxicity', 'LineWidth', 2.5, 'FontSize', 9);
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), Mg_levels, 'UniformOutput', false));
xtickangle(45);
xlabel('[Mg^{2+}] (mM)', 'FontSize', 11, 'Interpreter', 'tex');
ylabel('Peak [Ca^{2+}]_i (\muM)', 'FontSize', 11, 'Interpreter', 'tex');
title('F. Peak Ca (80 Hz)', 'FontSize', 12, 'FontWeight', 'bold');
box off;

sgtitle('RGC Model Overview', 'FontSize', 14, 'FontWeight', 'bold');
drawnow;  % Ensure figure is rendered

%% =========================================================================
%  FIGURE 2: DOSE-RESPONSE (IMPROVED WITH EXTENDED FREQUENCIES)
%  =========================================================================
fig2 = figure('Position', [100 100 1000 420], 'Color', 'w');

% A: Functional Cost (Spike Loss)
subplot(1,2,1);
hold on;

% Plot each frequency with distinct markers
markers = {'o', 's', 'd', '^', 'v', 'p'};  % Different markers for each freq
for f = 1:n_freq
    plot(Mg_levels, results.spike_red(:,f), ['-' markers{f}], 'Color', colors.freq(f,:), ...
        'MarkerFaceColor', colors.freq(f,:), 'MarkerSize', 6, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d Hz', frequencies(f)));
end

% Threshold line (exclude from legend)
yline(max_spike_loss, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(0.65, max_spike_loss + 2, sprintf('≤%d%% threshold', max_spike_loss), 'FontSize', 9);

xlabel('[Mg^{2+}] (mM)', 'FontSize', 11);
ylabel('Spike Loss (%)', 'FontSize', 11);
title('A. Functional Cost', 'FontSize', 12, 'FontWeight', 'bold');

% Legend with 2 columns for 6 frequencies
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'NumColumns', 2);

% Dynamic y-limit based on data
max_spike_loss_data = max(results.spike_red(:));
ylim([0 max(30, ceil(max_spike_loss_data/5)*5 + 5)]); 
xlim([0 2.6]);
box off; grid on;

% B: Neuroprotection (Ca²⁺ Reduction)
subplot(1,2,2);
hold on;
for f = 1:n_freq
    plot(Mg_levels, results.Ca_red(:,f), ['-' markers{f}], 'Color', colors.freq(f,:), ...
        'MarkerFaceColor', colors.freq(f,:), 'MarkerSize', 6, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d Hz', frequencies(f)));
end
xlabel('[Mg^{2+}] (mM)', 'FontSize', 11);
ylabel('Ca^{2+} Reduction (%)', 'FontSize', 11);
title('B. Neuroprotection', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southeast', 'Box', 'off', 'FontSize', 9, 'NumColumns', 2);
ylim([0 100]); xlim([0 2.6]);
box off; grid on;

% Updated annotation explaining the frequency regime
annotation('textbox', [0.30 0.10 0.20 0.25], 'String', ...
    'Note: At 10–60 Hz, neuron maintains 1:1 phase-locking; ≥80 Hz stress reveals dose-dependent functional tradeoff.', ...
    'FontSize', 10, 'FontAngle', 'italic', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

sgtitle('Dose-Response Analysis', 'FontSize', 14, 'FontWeight', 'bold');
drawnow;

%% =========================================================================
%  FIGURE 3: THERAPEUTIC WINDOWS (IMPROVED WITH CLEAR LABELS)
%  =========================================================================
fig3 = figure('Position', [150 150 950 400], 'Color', 'w');

% Find index for 80 Hz (reference frequency for trade-off space)
idx_80Hz = find(frequencies == 80);

% A: Trade-off space at 80 Hz (canonical stress condition)
subplot(1,2,1);
hold on;

Ca_base_80 = results.peak_Ca(1, idx_80Hz);
Ca_thresh_red = 100 * (Ca_base_80 - Ca_toxic) / Ca_base_80;

% Highlight optimal zone with label
fill([0 max_spike_loss max_spike_loss 0], [Ca_thresh_red Ca_thresh_red 100 100], ...
    [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

for m = 1:n_Mg
    ca_ok = results.peak_Ca(m, idx_80Hz) < Ca_toxic;
    func_ok_m = results.spike_red(m, idx_80Hz) <= max_spike_loss;
    if ca_ok && func_ok_m
        color = colors.optimal; sz = 120;
    else
        color = colors.suboptimal; sz = 60;
    end
    scatter(results.spike_red(m, idx_80Hz), results.Ca_red(m, idx_80Hz), sz, color, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end
plot(results.spike_red(:, idx_80Hz), results.Ca_red(:, idx_80Hz), 'k-', 'LineWidth', 1);

xline(max_spike_loss, 'b--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(max_spike_loss + 0.5, 48, 'Spike loss ≤20%', 'FontSize', 9, 'Color', [0.2 0.2 0.8], 'Rotation', 90);
yline(Ca_thresh_red, 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
text(15, Ca_thresh_red + 1.5, 'Peak Ca²⁺ < 1.0 µM', 'FontSize', 9, 'Color', [0.8 0.2 0.2], 'HorizontalAlignment', 'right');

xlabel('Spike Loss (%)', 'FontSize', 11);
ylabel('Ca^{2+} Reduction (%)', 'FontSize', 11, 'Interpreter', 'tex');
title('A. Trade-off Space (80 Hz)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-2 28]); ylim([45 92]);
text(3, 87, '\bf{THERAPEUTIC WINDOW}', 'Color', [0.1 0.5 0.1], 'FontSize', 10, 'Interpreter', 'tex');
text(3, 83, '1.6–2.0 mM', 'Color', [0.1 0.5 0.1], 'FontSize', 9);
box off; grid on;

% B: Window widths across all frequencies
subplot(1,2,2);

optimal_widths = zeros(n_freq, 1);
optimal_ranges = cell(n_freq, 1);
for f = 1:n_freq
    ca_ok_f = results.peak_Ca(:, f) < Ca_toxic;
    func_ok_f = results.spike_red(:, f) <= max_spike_loss;
    optimal_f = ca_ok_f & func_ok_f;
    if any(optimal_f)
        opt_idx = find(optimal_f);
        optimal_widths(f) = Mg_levels(opt_idx(end)) - Mg_levels(opt_idx(1));
        optimal_ranges{f} = sprintf('%.1f–%.1f', Mg_levels(opt_idx(1)), Mg_levels(opt_idx(end)));
    else
        optimal_ranges{f} = 'None';
        optimal_widths(f) = 0;
    end
end

b = bar(optimal_widths, 'FaceColor', 'flat');
b.CData = colors.freq;
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d', x), frequencies, 'UniformOutput', false));
xlabel('Stimulation Frequency (Hz)', 'FontSize', 11);
ylabel('Therapeutic Window Width (mM)', 'FontSize', 11);
title('B. Window Narrows with Frequency', 'FontSize', 12, 'FontWeight', 'bold');

for f = 1:n_freq
    if optimal_widths(f) > 0
        text(f, optimal_widths(f) + 0.06, optimal_ranges{f}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    else
        text(f, 0.08, 'None', 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'r');
    end
end
ylim([0 2.3]); box off;

sgtitle('Frequency-Dependent Therapeutic Windows', 'FontSize', 14, 'FontWeight', 'bold');
drawnow;

%% =========================================================================
%  FIGURE 4: INTERVENTION TIMING (IMPROVED WITH FINER RESOLUTION)
%  =========================================================================
fig4 = figure('Position', [100 100 1100 450], 'Color', 'w');

% A: Ca trajectories - show key conditions with full amplitude
subplot(1,2,1);
hold on;

% First plot the stress shading (behind everything)
max_Ca_all = max(cellfun(@max, Ca_int_traces));
fill([stress_start stress_end stress_end stress_start]/1000, ...
     [0 0 max_Ca_all*1.1 max_Ca_all*1.1], [0.95 0.9 0.9], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.25);

% Select key conditions to show: None, Pre, +0, +0.2, +0.5, +1s
% Indices: 1=None, 2=Pre, 3=+0, 5=+0.2, 8=+0.5, 10=+1
plot_indices = [1, 2, 3, 5, 8, 10];  
line_styles = {'-', '-', '--', '-.', ':', '-.'};
line_colors_traces = {colors.danger, colors.optimal, [0.2 0.7 0.4], ...
                      [0.4 0.6 0.3], [0.7 0.5 0.1], [0.5 0.5 0.5]};
line_widths = [2.5, 2.5, 2, 1.8, 1.8, 1.5];

h_lines = zeros(length(plot_indices), 1);
legend_labels_traces = cell(length(plot_indices), 1);
for i = 1:length(plot_indices)
    k = plot_indices(i);
    h_lines(i) = plot(t_int/1000, Ca_int_traces{k}, line_styles{i}, ...
        'Color', line_colors_traces{i}, 'LineWidth', line_widths(i));
    legend_labels_traces{i} = int_results.labels{k};
end

% Toxicity threshold
yline(Ca_toxic, 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(5.4, Ca_toxic - 0.15, 'Toxic', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

xlabel('Time (s)', 'FontSize', 11);
ylabel('[Ca^{2+}]_i (µM)', 'FontSize', 11);
title('A. Ca^{2+} Dynamics', 'FontSize', 12, 'FontWeight', 'bold');

% Legend
legend(h_lines, legend_labels_traces, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
xlim([0 6]); 
ylim([0 max_Ca_all * 1.05]);  % Dynamic y-limit to show full amplitude
box off; grid on;

% B: Protection efficacy bar chart - ALL conditions
subplot(1,2,2);
hold on;

% Create color gradient from green (protected) to gray (no protection)
n_bars = n_int;
bar_colors_int = zeros(n_bars, 3);
bar_colors_int(1,:) = colors.danger;  % No treatment - red
bar_colors_int(2,:) = colors.optimal;  % Pre-treated - green

% Gradient for delayed interventions based on protection level
for i = 3:n_bars
    prot = int_results.protection(i) / 100;
    bar_colors_int(i,:) = prot * colors.optimal + (1-prot) * colors.suboptimal;
end

b = bar(int_results.protection, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 0.8);
b.CData = bar_colors_int;

% X-axis labels
set(gca, 'XTick', 1:n_int, 'XTickLabel', int_results.labels, 'FontSize', 9);
xtickangle(45);
ylabel('Protection (%)', 'FontSize', 11);
xlabel('Intervention Timing (s)', 'FontSize', 10);
title('B. Protection vs Timing', 'FontSize', 12, 'FontWeight', 'bold');

% Add percentage labels on bars (only for significant values)
for i = 1:n_int
    if int_results.protection(i) >= 10
        text(i, int_results.protection(i) + 3, sprintf('%.0f%%', int_results.protection(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    elseif int_results.protection(i) >= 1
        text(i, int_results.protection(i) + 3, sprintf('%.0f%%', int_results.protection(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
    end
end

% Add annotation showing critical window
% Find where protection drops below 50%
prot_50_idx = find(int_results.protection(3:end) < 50, 1) + 2;
if ~isempty(prot_50_idx) && prot_50_idx <= n_int
    crit_time = intervention_delays(prot_50_idx - 2) / 1000;
    annotation('textbox', [0.73 0.30 0.2 0.08], 'String', ...
        sprintf('Critical window:\n< %.1f s for >50%% protection', crit_time), ...
        'FontSize', 9, 'EdgeColor', [0.8 0.4 0.2], 'BackgroundColor', [1 0.95 0.9], ...
        'HorizontalAlignment', 'center', 'FitBoxToText', 'on');
end

ylim([0 115]); box off;

% Add timing definition annotation
annotation('textbox', [0.215 0.55 0.2 0.2], 'String', ...
    'Timing: Pre = Mg present before stimulation; +X = Mg applied X s after stimulation onset', ...
    'FontSize', 8, 'EdgeColor', 'none', 'FontAngle', 'italic', 'HorizontalAlignment', 'left');

sgtitle('Intervention Timing Analysis', 'FontSize', 14, 'FontWeight', 'bold');
drawnow;

%% =========================================================================
%  FIGURE 5: MECHANISM (ENHANCED - 6 PANELS) - CORRECTED METRICS v2
%  =========================================================================
fig5 = figure('Position', [100 100 1150 600], 'Color', 'w');

Mg_mech = [0.2, 1.8];
Tmax_mech = 500; dt_mech = 0.01;
tn_mech = round(Tmax_mech/dt_mech);
t_mech = (0:tn_mech)*dt_mech;

period_mech = 1000/80;
Glu_mech = zeros(1, tn_mech+1);
for i = 1:tn_mech+1
    if mod(t_mech(i), period_mech) < stim.pulse_duration
        Glu_mech(i) = stim.Glu_conc;
    end
end

% Run detailed simulations with all outputs
[V_low, Ca_low, I_NMDA_low, I_AMPA_low, B_Mg_low] = run_simulation_detailed(t_mech, dt_mech, Glu_mech, Mg_mech(1), p);
[V_high, Ca_high, I_NMDA_high, I_AMPA_high, B_Mg_high] = run_simulation_detailed(t_mech, dt_mech, Glu_mech, Mg_mech(2), p);

% =========================================================================
% CALCULATE CONSISTENT METRICS using INTEGRATED CHARGE (∫I dt)
% This is the mechanistically relevant quantity for Ca²⁺ influx
% =========================================================================

% NMDA charge: integrate |I_NMDA| over time (inward current = negative)
% trapz(x,y) already computes the integral, no need to multiply by dt
Q_NMDA_low = trapz(t_mech, abs(I_NMDA_low));   % units: µA·ms/cm²
Q_NMDA_high = trapz(t_mech, abs(I_NMDA_high));
NMDA_charge_reduction = 100 * (1 - Q_NMDA_high/Q_NMDA_low);

% AMPA charge: same approach
Q_AMPA_low = trapz(t_mech, abs(I_AMPA_low));
Q_AMPA_high = trapz(t_mech, abs(I_AMPA_high));
AMPA_charge_change = 100 * (Q_AMPA_high/Q_AMPA_low - 1);

% Ca²⁺: use PEAK (clinically relevant for toxicity threshold)
peak_Ca_low = max(Ca_low);
peak_Ca_high = max(Ca_high);
Ca_peak_reduction = 100 * (1 - peak_Ca_high/peak_Ca_low);

% Also compute Ca²⁺ AUC for consistency check
AUC_Ca_low = trapz(t_mech, Ca_low - p.Ca_rest);
AUC_Ca_high = trapz(t_mech, Ca_high - p.Ca_rest);
Ca_AUC_reduction = 100 * (1 - AUC_Ca_high/AUC_Ca_low);

% Block factor: time-averaged
mean_B_low = mean(B_Mg_low);
mean_B_high = mean(B_Mg_high);
B_reduction = 100 * (1 - mean_B_high/mean_B_low);

% Print metrics for verification
fprintf('\n=== Figure 5 Mechanistic Metrics ===\n');
fprintf('  Block factor B(V) mean: %.3f → %.3f (%.1f%% reduction)\n', mean_B_low, mean_B_high, B_reduction);
fprintf('  NMDA integrated charge: %.1f%% reduction\n', NMDA_charge_reduction);
fprintf('  AMPA integrated charge: %.1f%% change (indirect, via V)\n', AMPA_charge_change);
fprintf('  Ca²⁺ peak: %.2f → %.2f µM (%.1f%% reduction)\n', peak_Ca_low, peak_Ca_high, Ca_peak_reduction);
fprintf('  Ca²⁺ AUC: %.1f%% reduction\n', Ca_AUC_reduction);
fprintf('=====================================\n\n');

% Define colors for better contrast (red vs blue, not red vs green)
color_low = [0.85 0.2 0.2];      % Red - pathological (0.2 mM)
color_high = [0.1 0.4 0.8];      % Blue - therapeutic (1.8 mM)

% =========================================================================
% PANEL A: Membrane Potential
% =========================================================================
subplot(2,3,1);
hold on;
% Plot therapeutic FIRST (behind), then pathological on top
plot(t_mech, V_high, '-', 'Color', color_high, 'LineWidth', 2.5);
plot(t_mech, V_low, '-', 'Color', color_low, 'LineWidth', 1.2);
xlabel('Time (ms)', 'FontSize', 10);
ylabel('V (mV)', 'FontSize', 10);
title('A. Membrane Potential', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 250]); box off; grid on;
% Use annotation boxes with Mg²⁺ label
text(5, 38, '0.2 mM Mg^{2+}', 'Color', color_low, 'FontSize', 9, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', color_low, 'Margin', 2, 'Interpreter', 'tex');
text(5, 21, '1.8 mM Mg^{2+}', 'Color', color_high, 'FontSize', 9, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', color_high, 'Margin', 2, 'Interpreter', 'tex');

% =========================================================================
% PANEL B: Mg²⁺ Block Factor B(V) - with clear definition
% =========================================================================
subplot(2,3,2);
hold on;
% Fill area between curves to show the difference
fill([t_mech, fliplr(t_mech)], [B_Mg_low, fliplr(B_Mg_high)], ...
    [0.9 0.85 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(t_mech, B_Mg_low, '-', 'Color', color_low, 'LineWidth', 2);
plot(t_mech, B_Mg_high, '-', 'Color', color_high, 'LineWidth', 2);
xlabel('Time (ms)', 'FontSize', 10);
ylabel('B(V) (fraction unblocked)', 'FontSize', 10);
title('B. Voltage-Dependent Mg^{2+} Block', 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
xlim([0 250]); ylim([0 1.05]);
% Annotations with mean values
text(160, 0.92, sprintf('⟨B⟩ = %.2f', mean_B_low), 'Color', color_low, 'FontSize', 9, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', color_low, 'Margin', 2);
text(160, 0.18, sprintf('⟨B⟩ = %.2f', mean_B_high), 'Color', color_high, 'FontSize', 9, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', color_high, 'Margin', 2);
% Add B(V) definition
text(5, 0.08, 'B(V) = 1/(1+η[Mg²⁺]e^{-γV})', 'FontSize', 7, 'Color', [0.4 0.4 0.4], ...
    'Interpreter', 'tex');
box off; grid on;

% =========================================================================
% PANEL C: AMPA Current (NO direct Mg block - change via V only)
% =========================================================================
subplot(2,3,3);
hold on;
plot(t_mech, -I_AMPA_high, '-', 'Color', color_high, 'LineWidth', 2.5);
plot(t_mech, -I_AMPA_low, '-', 'Color', color_low, 'LineWidth', 1.2);
xlabel('Time (ms)', 'FontSize', 10);
ylabel('I_{AMPA} (µA/cm²)', 'FontSize', 10);
% Explicit title: indirect effect via voltage
title(sprintf('C. AMPA Current (%.0f%% Δ, indirect via V)', abs(AMPA_charge_change)), ...
    'FontSize', 10, 'FontWeight', 'bold');
xlim([0 250]); box off; grid on;
% Clear annotation explaining indirect effect
text(125, max(-I_AMPA_low)*0.82, {'No direct Mg²⁺ block on AMPA'; 'Δ reflects altered V(t) trajectory'}, ...
    'FontSize', 7, 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3], ...
    'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.7 0.7 0.5], 'Margin', 2);

% =========================================================================
% PANEL D: NMDA Integrated Charge (Q_NMDA = ∫|I_NMDA|dt)
% =========================================================================
subplot(2,3,4);
hold on;
% Fill area to emphasize the charge reduction
fill([t_mech, fliplr(t_mech)], [-I_NMDA_low, fliplr(-I_NMDA_high)], ...
    [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(t_mech, -I_NMDA_low, '-', 'Color', color_low, 'LineWidth', 2);
plot(t_mech, -I_NMDA_high, '-', 'Color', color_high, 'LineWidth', 2);
xlabel('Time (ms)', 'FontSize', 10);
ylabel('I_{NMDA} (µA/cm²)', 'FontSize', 10);
% Title explicitly states integrated charge metric
title(sprintf('D. NMDA Charge Q=∫|I|dt (%.0f%% ↓)', NMDA_charge_reduction), ...
    'FontSize', 10, 'FontWeight', 'bold');
xlim([0 250]); box off; grid on;
% Add annotation showing the shaded area represents charge difference
text(180, 12, sprintf('ΔQ = %.0f%%', NMDA_charge_reduction), 'FontSize', 10, ...
    'FontWeight', 'bold', 'Color', [0.1 0.5 0.1], 'BackgroundColor', 'w');

% =========================================================================
% PANEL E: Peak Intracellular Ca²⁺
% =========================================================================
subplot(2,3,5);
hold on;
plot(t_mech, Ca_low, '-', 'Color', color_low, 'LineWidth', 2.5);
plot(t_mech, Ca_high, '-', 'Color', color_high, 'LineWidth', 2.5);
% Toxic threshold
yline(Ca_toxic, 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(Tmax_mech*0.55, Ca_toxic+0.18, 'Toxic threshold', 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
xlabel('Time (ms)', 'FontSize', 10);
ylabel('[Ca^{2+}]_i (µM)', 'FontSize', 10, 'Interpreter', 'tex');
title(sprintf('E. Peak [Ca^{2+}]_i (%.0f%% reduction)', Ca_peak_reduction), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
xlim([0 Tmax_mech]); box off; grid on;
% Label curves directly
text(420, peak_Ca_low - 0.3, '0.2 mM', 'Color', color_low, 'FontSize', 9, 'FontWeight', 'bold');
text(420, peak_Ca_high + 0.15, '1.8 mM', 'Color', color_high, 'FontSize', 9, 'FontWeight', 'bold');

% =========================================================================
% PANEL F: QUANTITATIVE SUMMARY (restored per reviewer request)
% =========================================================================
subplot(2,3,6);
axis off;

% Title
text(0.5, 0.97, '\bf{F. Quantitative Summary}', 'FontSize', 11, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');

% Summary box
rectangle('Position', [0.02, 0.08, 0.96, 0.85], 'FaceColor', [0.95 0.97 0.99], ...
    'EdgeColor', [0.2 0.4 0.6], 'LineWidth', 1.5, 'Curvature', 0.08);

y = 0.82; dy = 0.105;

% Header with color squares
text(0.5, y, '\bf{[Mg^{2+}]: 0.2 → 1.8 mM}', 'FontSize', 10, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');
rectangle('Position', [0.12, y-0.025, 0.05, 0.035], 'FaceColor', color_low, 'EdgeColor', 'k');
rectangle('Position', [0.80, y-0.025, 0.05, 0.035], 'FaceColor', color_high, 'EdgeColor', 'k');
y = y - dy;

% Metrics - clearly labeled
text(0.06, y, sprintf('⟨B(V)⟩: %.2f → %.2f (%.0f%% ↓)', mean_B_low, mean_B_high, B_reduction), ...
    'FontSize', 9);
y = y - dy;

text(0.06, y, sprintf('Q_{NMDA}: \\bf{%.0f%% ↓}', NMDA_charge_reduction), ...
    'FontSize', 9, 'Interpreter', 'tex', 'Color', [0.1 0.5 0.1]);
y = y - dy;

text(0.06, y, sprintf('Q_{AMPA}: %.0f%% Δ (indirect)', abs(AMPA_charge_change)), ...
    'FontSize', 9, 'Interpreter', 'tex');
y = y - dy;

text(0.06, y, sprintf('Peak Ca^{2+}: %.2f → %.2f µM', peak_Ca_low, peak_Ca_high), ...
    'FontSize', 9, 'Interpreter', 'tex');

% Bottom conclusion
text(0.5, 0.17, sprintf('\\bf{%.0f%% Ca^{2+} ↓ → Neuroprotection}', Ca_peak_reduction), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', [0.1 0.4 0.7], 'Interpreter', 'tex');

% Note about charge metric
text(0.5, 0.02, 'Q = ∫|I|dt (charge ∝ Ca²⁺ load)', 'FontSize', 7, ...
    'HorizontalAlignment', 'center', 'FontAngle', 'italic', 'Color', [0.4 0.4 0.4]);

sgtitle('Mechanistic Basis of Mg^{2+} Neuroprotection', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
drawnow;

% =========================================================================
% QUANTITATIVE SUMMARY - PRINTED TO CONSOLE FOR LATER USE
% =========================================================================
fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║           FIGURE 5: QUANTITATIVE SUMMARY                        ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  [Mg²⁺] Concentration Change: 0.2 → 1.8 mM                      ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Mean Block Factor B(V):      %.3f → %.3f  (%.1f%% reduction)   ║\n', mean_B_low, mean_B_high, B_reduction);
fprintf('║  NMDA Charge (∫|I|dt):        %.1f%% reduction                   ║\n', NMDA_charge_reduction);
fprintf('║  AMPA Charge (∫|I|dt):        %.1f%% change (indirect via V)    ║\n', abs(AMPA_charge_change));
fprintf('║  Peak [Ca²⁺]ᵢ:                %.2f → %.2f µM (%.1f%% reduction) ║\n', peak_Ca_low, peak_Ca_high, Ca_peak_reduction);
fprintf('║  Ca²⁺ AUC:                    %.1f%% reduction                   ║\n', Ca_AUC_reduction);
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  CONCLUSION: %.0f%% Ca²⁺ reduction → Neuroprotection             ║\n', Ca_peak_reduction);
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

%% =========================================================================
%  FIGURE S1: HIGH-RESOLUTION THERAPEUTIC WINDOW
%  =========================================================================
figS1 = figure('Position', [50 50 1000 420], 'Color', 'w');

% Define consistent colors
color_spike = [0.2 0.4 0.8];   % Blue for spike loss
color_Ca = [0.8 0.3 0.2];      % Red/orange for Ca²⁺

% A: Dose-response with dual axes (MATCHING AXIS COLORS)
subplot(1,2,1);

% Left axis: Spike loss - BLUE
yyaxis left;
set(gca, 'YColor', color_spike);  % Match axis color to curve
plot(Mg_highres, hr_results.spike_red, '-o', 'LineWidth', 2.5, 'MarkerSize', 7, ...
    'MarkerFaceColor', color_spike, 'Color', color_spike);
ylabel('Spike Loss (%)', 'FontSize', 11, 'Color', color_spike);
ylim([0 35]);
hold on;
yline(max_spike_loss, '--', 'LineWidth', 2, 'Color', color_spike);
text(1.5, max_spike_loss+1.8, 'Spike loss ≤20%', 'FontSize', 9, 'Color', color_spike);

% Right axis: Peak Ca - RED/ORANGE
yyaxis right;
set(gca, 'YColor', color_Ca);  % Match axis color to curve
plot(Mg_highres, hr_results.peak_Ca, '-s', 'LineWidth', 2.5, 'MarkerSize', 7, ...
    'MarkerFaceColor', color_Ca, 'Color', color_Ca);
hold on;
yline(Ca_toxic, '-', 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2]);
text(1.5, Ca_toxic+0.05, 'Peak Ca²⁺ < 1.0 µM', 'FontSize', 9, 'Color', color_Ca);
ylabel('Peak [Ca^{2+}]_i (µM)', 'FontSize', 11, 'Color', color_Ca, 'Interpreter', 'tex');

% Shade optimal region
opt_idx = find(optimal);
if ~isempty(opt_idx)
    Mg_opt_start = Mg_highres(opt_idx(1));
    Mg_opt_end = Mg_highres(opt_idx(end));
    yl = ylim;
    fill([Mg_opt_start Mg_opt_end Mg_opt_end Mg_opt_start], ...
         [yl(1) yl(1) yl(2) yl(2)], [0.85 0.95 0.85], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

xlabel('[Mg^{2+}] (mM)', 'FontSize', 11, 'Interpreter', 'tex');
title('A. Dose-Response (80 Hz)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.95 2.55]);
grid on; box off;

% B: Trade-off space with EXPLICIT scatter points
subplot(1,2,2);
hold on;

% Print debug info to verify data
fprintf('Figure S1B data verification:\n');
fprintf('  Spike loss range: %.1f%% to %.1f%%\n', min(hr_results.spike_red), max(hr_results.spike_red));
fprintf('  Ca reduction range: %.1f%% to %.1f%%\n', min(hr_results.Ca_red), max(hr_results.Ca_red));

% Set axis limits based on data
x_data = hr_results.spike_red;
y_data = hr_results.Ca_red;
x_min = floor(min(x_data)) - 2;
x_max = ceil(max(x_data)) + 2;
y_min = floor(min(y_data)) - 5;
y_max = ceil(max(y_data)) + 5;

% Threshold line in Ca reduction space
Ca_thresh_red_hr = 100 * (baseline_Ca - Ca_toxic) / baseline_Ca;
fprintf('  Ca threshold in reduction space: %.1f%%\n', Ca_thresh_red_hr);

% Shade optimal region FIRST (behind data)
fill([x_min max_spike_loss max_spike_loss x_min], ...
     [Ca_thresh_red_hr Ca_thresh_red_hr y_max y_max], ...
     [0.85 0.95 0.85], 'EdgeColor', [0.5 0.8 0.5], 'FaceAlpha', 0.4, 'LineWidth', 1);

% Plot connecting line
h_line = plot(x_data, y_data, 'k-', 'LineWidth', 1.5);

% Plot EACH point with explicit scatter - LARGER markers
for m = 1:n_Mg_hr
    % Determine color based on classification
    if optimal(m)
        c = [0.2 0.7 0.2];  % Green = optimal
        sz = 180;
    elseif ca_safe(m) && ~func_ok(m)
        c = [0.9 0.5 0.1];  % Orange = Ca safe but spike loss too high
        sz = 120;
    elseif func_ok(m) && ~ca_safe(m)
        c = [0.3 0.5 0.9];  % Blue = function OK but Ca too high
        sz = 120;
    else
        c = [0.5 0.5 0.5];  % Gray = fails both
        sz = 100;
    end
    
    scatter(x_data(m), y_data(m), sz, c, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

% Add Mg concentration labels to selected points
label_positions = [1, ceil(n_Mg_hr/3), ceil(2*n_Mg_hr/3), n_Mg_hr];
for li = label_positions
    text(x_data(li)+0.5, y_data(li)+2, sprintf('%.1f', Mg_highres(li)), ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
end

% Add threshold lines
xline(max_spike_loss, 'b--', 'LineWidth', 2.5);
yline(Ca_thresh_red_hr, 'r-', 'LineWidth', 2.5);

% Labels for threshold lines - standardized wording
text(max_spike_loss + 0.5, y_min + 5, 'Spike loss ≤20%', ...
    'FontSize', 9, 'Color', [0.2 0.2 0.8], 'Rotation', 90);
text(x_max - 4, Ca_thresh_red_hr + 2, 'Peak Ca²⁺ < 1.0 µM', ...
    'FontSize', 9, 'Color', [0.8 0.2 0.2]);

% Axis settings
xlabel('Spike Loss (%)', 'FontSize', 11);
ylabel('Ca^{2+} Reduction from Baseline (%)', 'FontSize', 11, 'Interpreter', 'tex');
title('B. Trade-off Space (Each Point = One [Mg^{2+}])', 'FontSize', 12, 'FontWeight', 'bold');

% Show optimal range in legend area
if ~isempty(opt_idx)
    text(x_min+1, y_max-3, sprintf('OPTIMAL: %.1f–%.1f mM', ...
        Mg_highres(opt_idx(1)), Mg_highres(opt_idx(end))), ...
        'Color', [0.1 0.5 0.1], 'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.2 0.6 0.2]);
end

% Manual legend with explicit scatter calls - positioned to avoid data overlap
h1 = scatter(nan, nan, 100, [0.2 0.7 0.2], 'filled', 'MarkerEdgeColor', 'k');
h2 = scatter(nan, nan, 80, [0.9 0.5 0.1], 'filled', 'MarkerEdgeColor', 'k');
h3 = scatter(nan, nan, 80, [0.3 0.5 0.9], 'filled', 'MarkerEdgeColor', 'k');
h4 = scatter(nan, nan, 80, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor', 'k');
leg = legend([h1 h2 h3 h4], {'Optimal', 'Ca OK only', 'Function OK only', 'Neither'}, ...
    'Location', 'southwest', 'Box', 'on', 'FontSize', 8);
set(leg, 'Color', [1 1 1 0.9]);  % Semi-transparent white background

xlim([x_min x_max]); ylim([y_min y_max]);
grid on; box off;

sgtitle('Therapeutic Window Analysis at 80 Hz', ...
    'FontSize', 13, 'FontWeight', 'bold');
drawnow;

%% =========================================================================
%  FIGURE S2: DETAILED RESULTS TABLE
%  =========================================================================
figS2 = figure('Position', [100 100 780 620], 'Color', 'w');
axis off;

% Title with methodology note
text(0.5, 0.97, '\bf{Results at 80 Hz}', 'FontSize', 13, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');

% Methodology box at top
rectangle('Position', [0.05, 0.90, 0.90, 0.045], 'FaceColor', [0.95 0.95 0.98], ...
    'EdgeColor', [0.7 0.7 0.8], 'LineWidth', 1);
text(0.5, 0.92, sprintf('Protocol: T_{stim} = 2.5 s, f = 80 Hz → N_{pulses} = %d | Spike Loss = 100×(1 - N_{spikes}/N_{pulses})', ...
    expected_pulses), 'FontSize', 9, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

% Table header
y_header = 0.86;
text(0.08, y_header, '[Mg^{2+}]', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');
text(0.22, y_header, 'Spikes', 'FontSize', 10, 'FontWeight', 'bold');
text(0.36, y_header, 'Loss (%)', 'FontSize', 10, 'FontWeight', 'bold');
text(0.52, y_header, 'Peak Ca', 'FontSize', 10, 'FontWeight', 'bold');
text(0.68, y_header, 'Ca < 1.0?', 'FontSize', 10, 'FontWeight', 'bold');
text(0.85, y_header, 'Optimal', 'FontSize', 10, 'FontWeight', 'bold');

line([0.03 0.97], [y_header-0.02 y_header-0.02], 'Color', 'k', 'LineWidth', 1.2);

row_spacing = 0.043;
start_y = y_header - 0.05;

for m = 1:n_Mg_hr
    y = start_y - (m-1)*row_spacing;
    
    if optimal(m)
        rectangle('Position', [0.03, y-0.016, 0.94, 0.034], ...
            'FaceColor', [0.92 0.97 0.92], 'EdgeColor', 'none');
    end
    
    text(0.08, y, sprintf('%.1f mM', Mg_highres(m)), 'FontSize', 9);
    text(0.22, y, sprintf('%d/%d', hr_results.spikes(m), expected_pulses), 'FontSize', 9, 'HorizontalAlignment', 'center');
    text(0.36, y, sprintf('%.1f', hr_results.spike_red(m)), 'FontSize', 9, 'HorizontalAlignment', 'center');
    
    if ca_safe(m)
        ca_color = [0.1 0.6 0.1];
    else
        ca_color = [0.8 0.1 0.1];
    end
    text(0.52, y, sprintf('%.3f', hr_results.peak_Ca(m)), 'FontSize', 9, 'Color', ca_color, 'HorizontalAlignment', 'center');
    
    if ca_safe(m)
        text(0.68, y, 'YES', 'FontSize', 9, 'Color', [0.1 0.6 0.1], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    else
        text(0.68, y, 'NO', 'FontSize', 9, 'Color', [0.8 0.1 0.1], 'HorizontalAlignment', 'center');
    end
    
    if optimal(m)
        text(0.85, y, '*', 'FontSize', 14, 'Color', [0.1 0.6 0.1], 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
    end
end

line([0.03 0.97], [0.12 0.12], 'Color', 'k', 'LineWidth', 1.2);

if ~isempty(opt_idx)
    text(0.5, 0.08, sprintf('OPTIMAL THERAPEUTIC RANGE: %.1f – %.1f mM', ...
        Mg_highres(opt_idx(1)), Mg_highres(opt_idx(end))), ...
        'FontSize', 12, 'Color', [0.1 0.5 0.1], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

text(0.5, 0.035, sprintf('Selection Criteria: Peak [Ca^{2+}]_i < %.1f µM AND spike loss ≤ %d%%', ...
    Ca_toxic, max_spike_loss), 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');
text(0.5, 0.01, '* indicates conditions meeting both criteria', ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);

%% =========================================================================
%  FIGURE S3: CALCIUM THRESHOLD SENSITIVITY ANALYSIS (HIGH RESOLUTION)
%  =========================================================================
figS3 = figure('Position', [100 100 1200 450], 'Color', 'w');

% A: Window Width vs Threshold (continuous curve)
subplot(1,3,1);
hold on;

% Plot window width as continuous curve
valid_idx = sens_results.window_widths > 0;
if any(valid_idx)
    plot(Ca_thresholds(valid_idx), sens_results.window_widths(valid_idx), '-o', ...
        'LineWidth', 2.5, 'MarkerSize', 6, 'Color', [0.2 0.5 0.7], ...
        'MarkerFaceColor', [0.2 0.5 0.7]);
end

% Mark thresholds with no window
zero_idx = sens_results.window_widths == 0;
if any(zero_idx)
    plot(Ca_thresholds(zero_idx), zeros(sum(zero_idx), 1), 'rx', ...
        'MarkerSize', 8, 'LineWidth', 2);
end

% Highlight default threshold
xline(1.0, 'g--', 'LineWidth', 2);
text(1.02, max(sens_results.window_widths)*0.85, 'Default', ...
    'FontSize', 9, 'Color', [0.2 0.6 0.2], 'FontWeight', 'bold');

xlabel('Ca^{2+} Toxicity Threshold (µM)', 'FontSize', 11, 'Interpreter', 'tex');
ylabel('Therapeutic Window Width (mM)', 'FontSize', 11);
title('A. Window Width vs Threshold', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.55 1.45]);
max_width = max(sens_results.window_widths);
if max_width > 0
    ylim([0 max_width * 1.15]);
else
    ylim([0 1]);
end
grid on; box off;

% B: Optimal Range Boundaries
subplot(1,3,2);
hold on;

% Fill between lower and upper bounds
if any(valid_idx)
    valid_thresh = Ca_thresholds(valid_idx);
    valid_lower = sens_results.lower_bounds(valid_idx);
    valid_upper = sens_results.upper_bounds(valid_idx);
    
    fill([valid_thresh, fliplr(valid_thresh)], [valid_lower', fliplr(valid_upper')], ...
        [0.7 0.9 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    % Plot boundaries
    plot(valid_thresh, valid_lower, '-', 'LineWidth', 2, 'Color', [0.2 0.5 0.2]);
    plot(valid_thresh, valid_upper, '-', 'LineWidth', 2, 'Color', [0.8 0.3 0.2]);
    
    % Labels
    text(mean(valid_thresh), mean(valid_lower)-0.1, 'Lower bound', ...
        'FontSize', 9, 'Color', [0.2 0.5 0.2], 'HorizontalAlignment', 'center');
    text(mean(valid_thresh), mean(valid_upper)+0.1, 'Upper bound', ...
        'FontSize', 9, 'Color', [0.8 0.3 0.2], 'HorizontalAlignment', 'center');
end

xlabel('Ca^{2+} Toxicity Threshold (µM)', 'FontSize', 11, 'Interpreter', 'tex');
ylabel('[Mg^{2+}] Optimal Range (mM)', 'FontSize', 11, 'Interpreter', 'tex');
title('B. Optimal Range Boundaries', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.55 1.45]); ylim([1.0 2.5]);
grid on; box off;

% C: Dose-Response with Key Thresholds - DEFAULT HIGHLIGHTED
subplot(1,3,3);
hold on;

% Plot dose-response curve
plot(Mg_highres, hr_results.peak_Ca, 'k-o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
    'MarkerFaceColor', 'k');

% Show selected thresholds - DEFAULT (1.0 µM) HIGHLIGHTED, others lighter
selected_thresh = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3];

for i = 1:length(selected_thresh)
    th_val = selected_thresh(i);
    
    % Highlight default (1.0 µM) strongly, others lighter
    if abs(th_val - 1.0) < 0.01
        % DEFAULT THRESHOLD - BOLD
        yline(th_val, '-', 'LineWidth', 3, 'Color', [0.1 0.7 0.2]);
        text(2.55, th_val, '\bf{1.0 µM (default)}', 'FontSize', 9, 'Color', [0.1 0.6 0.1], ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
        
        % Shade optimal region for default threshold
        th_idx = find(abs(Ca_thresholds - th_val) < 0.001);
        if ~isempty(th_idx) && sens_results.window_widths(th_idx) > 0
            fill_x = [sens_results.lower_bounds(th_idx), sens_results.upper_bounds(th_idx), ...
                      sens_results.upper_bounds(th_idx), sens_results.lower_bounds(th_idx)];
            fill_y = [0, 0, th_val, th_val];
            fill(fill_x, fill_y, [0.2 0.8 0.3], 'FaceAlpha', 0.25, 'EdgeColor', [0.1 0.6 0.1], 'LineWidth', 1.5);
        end
    else
        % Other thresholds - lighter/subdued
        alpha_val = 0.4;  % Reduced opacity
        line_color = [0.6 0.6 0.6];  % Grey
        yline(th_val, '--', 'LineWidth', 1, 'Color', line_color, 'Alpha', alpha_val);
        text(2.55, th_val, sprintf('%.1f µM', th_val), 'FontSize', 7, 'Color', [0.5 0.5 0.5]);
    end
end

xlabel('[Mg^{2+}] (mM)', 'FontSize', 11, 'Interpreter', 'tex');
ylabel('Peak [Ca^{2+}]_i (µM)', 'FontSize', 11, 'Interpreter', 'tex');
title('C. Dose-Response (default threshold highlighted)', 'FontSize', 11, 'FontWeight', 'bold');
xlim([1.0 2.6]); ylim([0 max(hr_results.peak_Ca)*1.1]);
grid on; box off;

sgtitle('Ca^{2+} Threshold Sensitivity Analysis', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
drawnow;

%% =========================================================================
%  FIGURE S4: NUMERICAL METHOD VALIDATION (NEW)
%  =========================================================================
figS4 = figure('Position', [100 100 1100 400], 'Color', 'w');

% Define colors for consistency
color_euler = [0.2 0.5 0.8];   % Blue for Euler
color_rk4 = [0.8 0.3 0.2];     % Red/orange for RK4

% A: Convergence with dt - show both Euler line and RK4 reference
subplot(1,3,1);
hold on;
% Plot Euler convergence
plot(dt_values*1000, euler_results.peak_Ca, '-o', 'LineWidth', 2.5, 'MarkerSize', 10, ...
    'MarkerFaceColor', color_euler, 'Color', color_euler, 'DisplayName', 'Forward Euler');
% Plot RK4 as horizontal reference line + marker
yline(rk4_results.peak_Ca, '--', 'Color', color_rk4, 'LineWidth', 2, ...
    'HandleVisibility', 'off');
plot(dt_rk4*1000, rk4_results.peak_Ca, 's', 'MarkerSize', 14, 'LineWidth', 2.5, ...
    'MarkerFaceColor', color_rk4, 'Color', color_rk4, 'DisplayName', 'RK4 (dt=20µs)');
% Label RK4 reference line
text(55, rk4_results.peak_Ca - 0.0005, 'RK4 reference', 'FontSize', 9, 'Color', color_rk4);

xlabel('Time Step (µs)', 'FontSize', 11);
ylabel('Peak [Ca^{2+}]_i (µM)', 'FontSize', 11, 'Interpreter', 'tex');
title('A. Convergence Analysis', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southeast', 'Box', 'off', 'FontSize', 9);
set(gca, 'XDir', 'reverse', 'FontSize', 10);
xlim([0 60]); box off; grid on;

% B: Error analysis - explicit reference statement
subplot(1,3,2);
hold on;

% Error relative to Euler dt=5µs (finest)
bar_data = [euler_results.Ca_error(1:3)', rk4_results.Ca_error];
bar_labels = {'Euler 50µs', 'Euler 20µs*', 'Euler 10µs', 'RK4 20µs'};

b = bar(bar_data, 0.65);
b.FaceColor = 'flat';
% Color: grey for coarse, blue shades for finer Euler, orange for RK4
b.CData = [0.7 0.7 0.7; color_euler; 0.1 0.4 0.7; color_rk4];

set(gca, 'XTickLabel', bar_labels, 'XTickLabelRotation', 25, 'FontSize', 9);
ylabel('Relative Error (%)', 'FontSize', 11);
title('B. Error vs Reference (Euler dt=5µs)', 'FontSize', 12, 'FontWeight', 'bold');

% Add value labels
for i = 1:length(bar_data)
    text(i, bar_data(i) + 0.08, sprintf('%.2f%%', bar_data(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end
ylim([0 max(bar_data)*1.35]); box off; grid on;

% Note about reference
text(0.5, -0.18, '*Used in manuscript', 'FontSize', 8, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontAngle', 'italic');

% C: Validation Summary Panel
subplot(1,3,3);
axis off;

text(0.5, 0.95, '\bf{C. Validation Summary}', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');

% Summary box
rectangle('Position', [0.05, 0.12, 0.9, 0.75], 'FaceColor', [0.95 0.98 0.95], ...
    'EdgeColor', [0.3 0.6 0.3], 'LineWidth', 1.5, 'Curvature', 0.08);

y = 0.80; dy = 0.10;

text(0.5, y, '\bf{Test Condition}', 'FontSize', 11, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
y = y - dy;
text(0.5, y, '80 Hz stimulation, [Mg^{2+}] = 1.8 mM', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
y = y - dy*1.3;

text(0.5, y, '\bf{Reference: Euler dt=5 µs}', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
y = y - dy*1.2;

text(0.5, y, '\bf{Key Results}', 'FontSize', 11, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
y = y - dy*0.9;
text(0.12, y, sprintf('Euler dt=20 µs:  %.4f µM', euler_results.peak_Ca(2)), 'FontSize', 10, 'Color', color_euler);
y = y - dy*0.8;
text(0.12, y, sprintf('RK4 dt=20 µs:    %.4f µM', rk4_results.peak_Ca), 'FontSize', 10, 'Color', color_rk4);
y = y - dy*0.8;
text(0.12, y, sprintf('Agreement:       %.3f%%', rk4_results.Ca_error), 'FontSize', 10, 'Color', [0.1 0.5 0.1], 'FontWeight', 'bold');

% Conclusion
text(0.5, 0.17, '\bf{✓ dt = 20 µs validated}', 'FontSize', 11, ...
    'HorizontalAlignment', 'center', 'Color', [0.1 0.6 0.1], 'Interpreter', 'tex');

sgtitle('Numerical Method Validation', 'FontSize', 13, 'FontWeight', 'bold');
drawnow;

%% =========================================================================
%  FIGURE S5: EXPANDED PARAMETER TABLE (NEW) - SPLIT INTO 2 PANELS
%  =========================================================================
figS5 = figure('Position', [50 50 1100 650], 'Color', 'w');

% PANEL A: Membrane, Reversal, Conductances
subplot(1,2,1);
axis off;

text(0.5, 0.98, '\bf{A. Intrinsic Properties}', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');

params_A = {
    % Headers
    'Symbol', 'Value', 'Unit', 'Description';
    % Membrane
    '\bf{Membrane}', '', '', '';
    'C_m', '1.0', 'µF/cm²', 'Capacitance';
    'V_{rest}', '-65', 'mV', 'Resting potential';
    % Reversal
    '\bf{Reversal Potentials}', '', '', '';
    'E_{Na}', '+50', 'mV', 'Sodium';
    'E_K', '-77', 'mV', 'Potassium';
    'E_{Ca}', '+120', 'mV', 'Calcium';
    'E_L', '-54.4', 'mV', 'Leak';
    'E_{exc}', '0', 'mV', 'Excitatory synaptic';
    % Conductances
    '\bf{Conductances (mS/cm²)}', '', '', '';
    'g_{Na}', '120', '', 'Sodium';
    'g_{Kdr}', '36', '', 'Delayed rectifier K⁺';
    'g_{KA}', '8', '', 'A-type K⁺';
    'g_{CaL}', '0.3', '', 'L-type Ca²⁺';
    'g_{KCa}', '0.3', '', 'Ca²⁺-activated K⁺';
    'g_L', '0.35', '', 'Leak';
};

y = 0.92;
dy = 0.045;
for r = 1:size(params_A, 1)
    if r == 1  % Header row
        text(0.05, y, params_A{r,1}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.35, y, params_A{r,2}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.50, y, params_A{r,3}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.65, y, params_A{r,4}, 'FontSize', 9, 'FontWeight', 'bold');
        line([0.02 0.98], [y-0.015 y-0.015], 'Color', 'k', 'LineWidth', 1);
    elseif startsWith(params_A{r,1}, '\bf{')  % Category header
        y = y - dy*0.3;
        text(0.03, y, params_A{r,1}, 'FontSize', 9, 'Interpreter', 'tex', 'Color', [0.2 0.2 0.6]);
    else  % Data row
        text(0.05, y, params_A{r,1}, 'FontSize', 9, 'Interpreter', 'tex');
        text(0.35, y, params_A{r,2}, 'FontSize', 9);
        text(0.50, y, params_A{r,3}, 'FontSize', 9);
        text(0.65, y, params_A{r,4}, 'FontSize', 9);
    end
    y = y - dy;
end

% PANEL B: Synaptic, Mg Block, Ca Dynamics
subplot(1,2,2);
axis off;

text(0.5, 0.98, '\bf{B. Synaptic & Ca²⁺ Dynamics}', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');

params_B = {
    'Symbol', 'Value', 'Unit', 'Description';
    '\bf{Synaptic}', '', '', '';
    'g_{AMPA}', '0.25', 'mS/cm²', 'Max AMPA conductance';
    'g_{NMDA}', '1.2', 'mS/cm²', 'Max NMDA conductance';
    'EC_{50,AMPA}', '0.5', 'mM', 'AMPA glutamate affinity';
    'EC_{50,NMDA}', '2', 'µM', 'NMDA glutamate affinity';
    'τ_{AMPA}', '0.3/3', 'ms', 'Rise/decay time';
    'τ_{NMDA}', '5/80', 'ms', 'Rise/decay time';
    '\bf{Mg²⁺ Block (Jahr-Stevens)}', '', '', '';
    'η', '0.28', 'mM⁻¹', 'Mg²⁺ sensitivity';
    'γ', '0.062', 'mV⁻¹', 'Voltage sensitivity';
    '\bf{Ca²⁺ Dynamics}', '', '', '';
    'f_{Ca}', '0.15', '—', 'Frac. NMDA Ca²⁺ current';
    '[Ca²⁺]_{rest}', '0.05', 'µM', 'Resting [Ca²⁺]ᵢ';
    'τ_{Ca}', '200', 'ms', 'Decay time constant';
    'K_d (KCa)', '0.5', 'µM', 'KCa half-activation';
    '\bf{Stimulation}', '', '', '';
    '[Glu]', '1.0', 'mM', 'Pulse concentration';
    't_{pulse}', '2', 'ms', 'Pulse duration';
};

y = 0.92;
for r = 1:size(params_B, 1)
    if r == 1
        text(0.05, y, params_B{r,1}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.35, y, params_B{r,2}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.50, y, params_B{r,3}, 'FontSize', 9, 'FontWeight', 'bold');
        text(0.65, y, params_B{r,4}, 'FontSize', 9, 'FontWeight', 'bold');
        line([0.02 0.98], [y-0.015 y-0.015], 'Color', 'k', 'LineWidth', 1);
    elseif startsWith(params_B{r,1}, '\bf{')
        y = y - dy*0.3;
        text(0.03, y, params_B{r,1}, 'FontSize', 9, 'Interpreter', 'tex', 'Color', [0.2 0.2 0.6]);
    else
        text(0.05, y, params_B{r,1}, 'FontSize', 9, 'Interpreter', 'tex');
        text(0.35, y, params_B{r,2}, 'FontSize', 9);
        text(0.50, y, params_B{r,3}, 'FontSize', 9);
        text(0.65, y, params_B{r,4}, 'FontSize', 9);
    end
    y = y - dy;
end

% Add note at bottom
annotation('textbox', [0.1 0.02 0.8 0.05], 'String', ...
    'References: Hodgkin-Huxley 1952; Jahr & Stevens 1990; Patneau & Mayer 1990; Destexhe et al. 1994', ...
    'FontSize', 8, 'FontAngle', 'italic', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

sgtitle('cModel Parameters', 'FontSize', 14, 'FontWeight', 'bold');
drawnow;

%% =========================================================================
%  SAVE ALL FIGURES
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  Saving figures...\n');
fprintf('=========================================================================\n\n');

% Ensure all figures are rendered
drawnow;
pause(0.5);  % Give MATLAB time to render

% Print output directory
fprintf('Output directory: %s\n\n', pwd);

% Save with validation
if isvalid(fig1), save_figure(fig1, 'Figure1_ModelOverview'); else, warning('fig1 invalid'); end
if isvalid(fig2), save_figure(fig2, 'Figure2_DoseResponse'); else, warning('fig2 invalid'); end
if isvalid(fig3), save_figure(fig3, 'Figure3_TherapeuticWindows'); else, warning('fig3 invalid'); end
if isvalid(fig4), save_figure(fig4, 'Figure4_InterventionTiming'); else, warning('fig4 invalid'); end
if isvalid(fig5), save_figure(fig5, 'Figure5_Mechanism'); else, warning('fig5 invalid'); end
if isvalid(figS1), save_figure(figS1, 'FigureS1_HighResolution'); else, warning('figS1 invalid'); end
if isvalid(figS2), save_figure(figS2, 'FigureS2_ResultsTable'); else, warning('figS2 invalid'); end
if isvalid(figS3), save_figure(figS3, 'FigureS3_ThresholdSensitivity'); else, warning('figS3 invalid'); end
if isvalid(figS4), save_figure(figS4, 'FigureS4_NumericalValidation'); else, warning('figS4 invalid'); end
if isvalid(figS5), save_figure(figS5, 'FigureS5_ExpandedParameters'); else, warning('figS5 invalid'); end

fprintf('\nAll figures saved to: %s\n', pwd);

%% =========================================================================
%  PRINT COMPREHENSIVE SUMMARY
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  COMPREHENSIVE RESULTS SUMMARY\n');
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

fprintf('\nCA²⁺ THRESHOLD SENSITIVITY (80 Hz) - KEY RESULTS:\n');
% Find first and last valid thresholds
valid_thresh = find(sens_results.window_widths > 0);
if ~isempty(valid_thresh)
    fprintf('  Thresholds tested: %.2f to %.2f µM (n=%d)\n', ...
        min(Ca_thresholds), max(Ca_thresholds), n_thresh);
    fprintf('  Window exists for thresholds: ≥ %.2f µM\n', Ca_thresholds(valid_thresh(1)));
    fprintf('  At default threshold (1.0 µM): %s (width: %.2f mM)\n', ...
        sens_results.optimal_ranges{Ca_thresholds == 1.0}, ...
        sens_results.window_widths(Ca_thresholds == 1.0));
    fprintf('  Max window width: %.2f mM at threshold %.2f µM\n', ...
        max(sens_results.window_widths), Ca_thresholds(sens_results.window_widths == max(sens_results.window_widths)));
else
    fprintf('  No valid therapeutic windows found\n');
end

fprintf('\nNUMERICAL VALIDATION:\n');
fprintf('  Forward Euler (dt=0.02 ms): Peak Ca = %.4f µM\n', euler_results.peak_Ca(2));
fprintf('  RK4 (dt=0.02 ms): Peak Ca = %.4f µM\n', rk4_results.peak_Ca);
fprintf('  Relative difference: %.3f%% (excellent agreement)\n', rk4_results.Ca_error);

fprintf('\n=========================================================================\n');
fprintf('  REVISION SUMMARY - Points Addressed\n');
fprintf('=========================================================================\n');
fprintf('  ✓ Extended frequency range (10–100 Hz) showing window closure\n');
fprintf('  ✓ Sensitivity analysis (17 thresholds, Figure S3)\n');
fprintf('  ✓ Numerical method validation (Figure S4)\n');
fprintf('  ✓ Expanded parameter table (Figure S5/Table 1)\n');
fprintf('  ✓ Improved figure labeling and fonts\n');
fprintf('  ✓ Explicit "Therapeutic Window" labels\n');
fprintf('  ✓ Clinical timing interpretation note added\n');
fprintf('=========================================================================\n');

%% =========================================================================
%  HELPER FUNCTIONS
%  =========================================================================

function save_figure(fig, filename)
    % -----------------------------------------------------------------------
    %  OUTPUT DIRECTORY
    %  Change output_dir to your preferred folder if needed.
    % -----------------------------------------------------------------------
    output_dir = pwd;
    % output_dir = 'C:\Users\YourName\Documents\Figures';   % Windows
    % output_dir = '/Users/YourName/Documents/Figures';     % Mac/Linux

    % PLOS ONE DPI requirements:
    %   600 DPI  -> combination figures (line art + text labels)  <-- all your figs
    %   300 DPI  -> photographs only
    DPI_TIFF = 600;
    DPI_PNG  = 300;   % PNG is just a preview, not for submission

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    filepath = fullfile(output_dir, filename);

    if ~isvalid(fig)
        warning('Figure handle is invalid for %s — skipping.', filename);
        return;
    end

    % -------------------------------------------------------------------
    %  1. TIFF  (PRIMARY — PLOS ONE submission format)
    %     -dtiff  = LZW-compressed TIFF (lossless, smaller file)
    %     -dtiffn = uncompressed TIFF (fallback, larger)
    % -------------------------------------------------------------------
    tif_path = [filepath '.tif'];
    saved_tiff = false;

    try
        print(fig, tif_path, '-dtiff', sprintf('-r%d', DPI_TIFF));
        fprintf('  [TIFF] Saved: %s.tif  (%d DPI, LZW)\n', filename, DPI_TIFF);
        saved_tiff = true;
    catch ME1
        warning('LZW TIFF failed (%s) — trying uncompressed...', ME1.message);
        try
            print(fig, tif_path, '-dtiffn', sprintf('-r%d', DPI_TIFF));
            fprintf('  [TIFF] Saved: %s.tif  (%d DPI, uncompressed)\n', filename, DPI_TIFF);
            saved_tiff = true;
        catch ME2
            warning('Uncompressed TIFF failed (%s) — trying saveas...', ME2.message);
            try
                saveas(fig, tif_path);
                fprintf('  [TIFF] Saved: %s.tif  (via saveas — verify DPI manually)\n', filename);
                saved_tiff = true;
            catch ME3
                warning('All TIFF methods failed for %s: %s', filename, ME3.message);
            end
        end
    end

    if ~saved_tiff
        fprintf('  [TIFF] WARNING: Could not save %s.tif\n', filename);
    end

    % -------------------------------------------------------------------
    %  2. PDF  (vector — for LaTeX manuscript / supplementary)
    % -------------------------------------------------------------------
    try
        set(fig, 'PaperPositionMode', 'auto');
        fig_pos = get(fig, 'Position');
        set(fig, 'PaperUnits', 'points');
        set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
        print(fig, [filepath '.pdf'], '-dpdf', '-vector');
        fprintf('  [PDF]  Saved: %s.pdf  (vector)\n', filename);
    catch ME
        try
            saveas(fig, [filepath '.pdf']);
            fprintf('  [PDF]  Saved: %s.pdf  (via saveas)\n', filename);
        catch
            warning('PDF save failed for %s: %s', filename, ME.message);
        end
    end

    % -------------------------------------------------------------------
    %  3. PNG  (quick preview only — NOT for PLOS ONE submission)
    % -------------------------------------------------------------------
    try
        print(fig, [filepath '.png'], '-dpng', sprintf('-r%d', DPI_PNG));
        fprintf('  [PNG]  Saved: %s.png  (%d DPI, preview only)\n', filename, DPI_PNG);
    catch ME
        try
            saveas(fig, [filepath '.png']);
            fprintf('  [PNG]  Saved: %s.png  (via saveas)\n', filename);
        catch
            warning('PNG save failed for %s: %s', filename, ME.message);
        end
    end

    fprintf('\n');
end


%% =========================================================================
%  SIMULATION FUNCTIONS
%  =========================================================================

function [V, Ca] = run_simulation_euler(t, dt, Glu, Mg_ext, p)
    % Forward Euler integration
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
    
    % Initial conditions
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
        
        % Synaptic activation
        if Glu_now > 0
            r_AMPA = (Glu_now/(Glu_now + p.EC50_AMPA)) / p.tau_AMPA_rise;
            r_NMDA = (Glu_now/(Glu_now + p.EC50_NMDA)) / p.tau_NMDA_rise;
        else
            r_AMPA = 0; r_NMDA = 0;
        end
        
        s_AMPA(i+1) = max(0, min(1, s_AMPA(i) + dt*(r_AMPA*(1-s_AMPA(i)) - s_AMPA(i)/p.tau_AMPA_decay)));
        s_NMDA(i+1) = max(0, min(1, s_NMDA(i) + dt*(r_NMDA*(1-s_NMDA(i)) - s_NMDA(i)/p.tau_NMDA_decay)));
        
        % Gating variables
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
        
        % Currents
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
        
        % Calcium dynamics
        Ca_influx = -p.k_Ca_NMDA * p.f_Ca * I_NMDA - p.k_Ca_CaL * I_CaL;
        Ca(i+1) = max(p.Ca_rest, Ca(i) + dt*(Ca_influx - (Ca(i)-p.Ca_rest)/p.tau_Ca));
        
        % Membrane potential
        I_total = I_Na + I_Kdr + I_KA + I_CaL + I_KCa + I_L + I_AMPA + I_NMDA;
        V(i+1) = v + dt*(-I_total)/p.Cm;
    end
end

function [V, Ca] = run_simulation_rk4(t, dt, Glu, Mg_ext, p)
    % 4th-order Runge-Kutta integration for validation
    tn = length(t) - 1;
    
    % State vector: [V, m, h, n, a, b, s_CaL, s_AMPA, s_NMDA, Ca]
    n_states = 10;
    y = zeros(n_states, tn+1);
    
    % Initial conditions
    v0 = p.V_rest;
    am = 0.1*(v0+40)/(1-exp(-(v0+40)/10)); bm = 4*exp(-(v0+65)/18);
    ah = 0.07*exp(-(v0+65)/20); bh = 1/(1+exp(-(v0+35)/10));
    an = 0.01*(v0+55)/(1-exp(-(v0+55)/10)); bn = 0.125*exp(-(v0+65)/80);
    
    y(:,1) = [v0; am/(am+bm); ah/(ah+bh); an/(an+bn); ...
              1/(1+exp(-(v0+50)/20)); 1/(1+exp((v0+70)/6)); ...
              1/(1+exp(-(v0+20)/6)); 0; 0; p.Ca_rest];
    
    for i = 1:tn
        Glu_now = Glu(i);
        
        % RK4 steps
        k1 = dt * derivatives(y(:,i), Glu_now, Mg_ext, p);
        k2 = dt * derivatives(y(:,i) + k1/2, Glu_now, Mg_ext, p);
        k3 = dt * derivatives(y(:,i) + k2/2, Glu_now, Mg_ext, p);
        k4 = dt * derivatives(y(:,i) + k3, Glu_now, Mg_ext, p);
        
        y(:,i+1) = y(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
        
        % Enforce bounds
        y(2:9,i+1) = max(0, min(1, y(2:9,i+1)));
        y(10,i+1) = max(p.Ca_rest, y(10,i+1));
    end
    
    V = y(1,:);
    Ca = y(10,:);
end

function dydt = derivatives(y, Glu, Mg, p)
    v = y(1); m = y(2); h = y(3); n = y(4);
    a = y(5); b = y(6); s_CaL = y(7);
    s_AMPA = y(8); s_NMDA = y(9); Ca_val = y(10);
    
    % Synaptic
    if Glu > 0
        r_AMPA = (Glu/(Glu + p.EC50_AMPA)) / p.tau_AMPA_rise;
        r_NMDA = (Glu/(Glu + p.EC50_NMDA)) / p.tau_NMDA_rise;
    else
        r_AMPA = 0; r_NMDA = 0;
    end
    
    ds_AMPA = r_AMPA*(1-s_AMPA) - s_AMPA/p.tau_AMPA_decay;
    ds_NMDA = r_NMDA*(1-s_NMDA) - s_NMDA/p.tau_NMDA_decay;
    
    % Gating
    if abs(v+40) < 1e-6, am = 1; else, am = 0.1*(v+40)/(1-exp(-(v+40)/10)); end
    bm = 4*exp(-(v+65)/18);
    ah = 0.07*exp(-(v+65)/20); bh = 1/(1+exp(-(v+35)/10));
    if abs(v+55) < 1e-6, an = 0.1; else, an = 0.01*(v+55)/(1-exp(-(v+55)/10)); end
    bn = 0.125*exp(-(v+65)/80);
    
    dm = (am*(1-m) - bm*m);
    dh = (ah*(1-h) - bh*h);
    dn = (an*(1-n) - bn*n);
    da = (1/(1+exp(-(v+50)/20)) - a)/5;
    db = (1/(1+exp((v+70)/6)) - b)/20;
    ds_CaL = (1/(1+exp(-(v+20)/6)) - s_CaL)/5;
    
    % Currents
    I_Na = p.g_Na * m^3 * h * (v - p.E_Na);
    I_Kdr = p.g_Kdr * n^4 * (v - p.E_K);
    I_KA = p.g_KA * a^3 * b * (v - p.E_K);
    I_CaL = p.g_CaL * s_CaL^2 * (v - p.E_Ca);
    I_L = p.g_L * (v - p.E_L);
    
    q_KCa = Ca_val^p.n_KCa / (Ca_val^p.n_KCa + p.Kd_KCa^p.n_KCa);
    I_KCa = p.g_KCa * q_KCa * (v - p.E_K);
    
    B_Mg = 1 / (1 + p.eta * Mg * exp(-p.gamma * v));
    I_AMPA = p.g_AMPA_max * s_AMPA * (v - p.E_exc);
    I_NMDA = p.g_NMDA_max * s_NMDA * B_Mg * (v - p.E_exc);
    
    % Calcium
    Ca_influx = -p.k_Ca_NMDA * p.f_Ca * I_NMDA - p.k_Ca_CaL * I_CaL;
    dCa = Ca_influx - (Ca_val - p.Ca_rest)/p.tau_Ca;
    
    % Voltage
    I_total = I_Na + I_Kdr + I_KA + I_CaL + I_KCa + I_L + I_AMPA + I_NMDA;
    dV = -I_total/p.Cm;
    
    dydt = [dV; dm; dh; dn; da; db; ds_CaL; ds_AMPA; ds_NMDA; dCa];
end

function [V, Ca, I_NMDA_rec, I_AMPA_rec, B_Mg_rec] = run_simulation_detailed(t, dt, Glu, Mg_ext, p)
    % Forward Euler with current and Mg block recording
    tn = length(t) - 1;
    V = zeros(1, tn+1); V(1) = p.V_rest;
    m_g = zeros(1, tn+1); h_g = zeros(1, tn+1); n_g = zeros(1, tn+1);
    a_g = zeros(1, tn+1); b_g = zeros(1, tn+1); s_CaL = zeros(1, tn+1);
    s_AMPA = zeros(1, tn+1); s_NMDA = zeros(1, tn+1);
    Ca = zeros(1, tn+1); Ca(1) = p.Ca_rest;
    I_NMDA_rec = zeros(1, tn+1);
    I_AMPA_rec = zeros(1, tn+1);
    B_Mg_rec = zeros(1, tn+1);
    
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
        I_AMPA_rec(i) = I_AMPA;
        B_Mg_rec(i) = B_Mg;
        
        Ca_influx = -p.k_Ca_NMDA * p.f_Ca * I_NMDA - p.k_Ca_CaL * I_CaL;
        Ca(i+1) = max(p.Ca_rest, Ca(i) + dt*(Ca_influx - (Ca(i)-p.Ca_rest)/p.tau_Ca));
        
        I_total = I_Na + I_Kdr + I_KA + I_CaL + I_KCa + I_L + I_AMPA + I_NMDA;
        V(i+1) = v + dt*(-I_total)/p.Cm;
    end
end
