%% =========================================================================
%  Mg²⁺ NEUROPROTECTION - SUPPLEMENTARY ANALYSES
%  =========================================================================
%
%  Figure S1: High-Resolution Therapeutic Window (80 Hz)
%  Figure S2: Detailed Results Tables
%
%  Author: Mehdi Borjkhani
%  Institution: ICTER, Polish Academy of Sciences
%  =========================================================================

clear; clc; close all;

fprintf('=========================================================================\n');
fprintf('  SUPPLEMENTARY ANALYSES\n');
fprintf('=========================================================================\n\n');

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
%  FIGURE S1: HIGH-RESOLUTION THERAPEUTIC WINDOW
%  =========================================================================
fprintf('Running high-resolution analysis...\n');

Mg_highres = 1.0:0.1:2.5;
n_Mg_hr = length(Mg_highres);

Tmax = 3000; dt = 0.02; transient = 500;
tn = round(Tmax/dt);
t = (0:tn)*dt;
idx_ss = t >= transient;

period = 1000/80;
expected_pulses = floor((Tmax - transient) / period);
Glu = zeros(1, tn+1);
for i = 1:tn+1
    if mod(t(i), period) < 2, Glu(i) = 1.0; end
end

hr_results.spikes = zeros(n_Mg_hr, 1);
hr_results.peak_Ca = zeros(n_Mg_hr, 1);

for m = 1:n_Mg_hr
    [V, Ca] = run_sim(t, dt, Glu, Mg_highres(m), p);
    V_ss = V(idx_ss); Ca_ss = Ca(idx_ss);
    hr_results.spikes(m) = sum(V_ss(1:end-1) < -20 & V_ss(2:end) >= -20);
    hr_results.peak_Ca(m) = max(Ca_ss);
end

baseline_spikes = 200;
hr_results.spike_red = 100 * (baseline_spikes - hr_results.spikes) / baseline_spikes;
hr_results.Ca_red = 100 * (4.59 - hr_results.peak_Ca) / 4.59;

ca_safe = hr_results.peak_Ca < Ca_toxic;
func_ok = hr_results.spike_red <= max_spike_loss;
optimal = ca_safe & func_ok;

fprintf('Done.\n\n');

% Create figure
figS1 = figure('Position', [50 50 900 350], 'Color', 'w');

% A: Dose-response
subplot(1,2,1);
yyaxis left;
plot(Mg_highres, hr_results.spike_red, '-o', 'LineWidth', 2, 'MarkerSize', 5);
ylabel('Spike Reduction (%)');
ylim([0 28]);

yyaxis right;
plot(Mg_highres, hr_results.peak_Ca, '-s', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
yline(Ca_toxic, 'r-', 'Toxic', 'LineWidth', 2);
ylabel('Peak Ca²⁺ (µM)');

yyaxis left;
yline(max_spike_loss, 'b--', sprintf('≤%d%%', max_spike_loss), 'LineWidth', 1.5);

xlabel('[Mg²⁺] (mM)');
title('A. High-Resolution Dose-Response');
grid on; box off;

% B: Optimal zone
subplot(1,2,2);
hold on;

colors_hr = zeros(n_Mg_hr, 3);
for m = 1:n_Mg_hr
    if optimal(m)
        colors_hr(m,:) = [0.2 0.7 0.2];
    elseif ca_safe(m)
        colors_hr(m,:) = [0.9 0.6 0.1];
    elseif func_ok(m)
        colors_hr(m,:) = [0.3 0.5 0.8];
    else
        colors_hr(m,:) = [0.6 0.6 0.6];
    end
end

Ca_thresh_red = 100 * (4.59 - Ca_toxic) / 4.59;

fill([0 max_spike_loss max_spike_loss 0], [Ca_thresh_red Ca_thresh_red 100 100], ...
    [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

for m = 1:n_Mg_hr
    sz = 60;
    if optimal(m), sz = 100; end
    scatter(hr_results.spike_red(m), hr_results.Ca_red(m), sz, colors_hr(m,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

plot(hr_results.spike_red, hr_results.Ca_red, 'k-', 'LineWidth', 0.5);

xline(max_spike_loss, 'b--', 'LineWidth', 1.5);
yline(Ca_thresh_red, 'r-', 'LineWidth', 2);

xlabel('Spike Reduction (%)');
ylabel('Ca²⁺ Reduction (%)');
title('B. Optimal Zone');

opt_idx = find(optimal);
if ~isempty(opt_idx)
    text(12, 85, sprintf('Optimal: %.1f–%.1f mM', Mg_highres(opt_idx(1)), Mg_highres(opt_idx(end))), ...
        'FontWeight', 'bold', 'Color', [0.2 0.5 0.2], 'FontSize', 9);
end

xlim([8 26]); ylim([68 88]);
grid on; box off;

sgtitle('Figure S1: High-Resolution Analysis (80 Hz)', 'FontSize', 12, 'FontWeight', 'bold');

%% =========================================================================
%  FIGURE S2: RESULTS TABLE
%  =========================================================================
figS2 = figure('Position', [100 100 700 550], 'Color', 'w');

axis off;

% Title
text(0.5, 0.95, 'High-Resolution Results at 80 Hz', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

% Table header
text(0.10, 0.88, '[Mg²⁺]', 'FontSize', 9, 'FontWeight', 'bold');
text(0.25, 0.88, 'Spikes', 'FontSize', 9, 'FontWeight', 'bold');
text(0.40, 0.88, 'Loss (%)', 'FontSize', 9, 'FontWeight', 'bold');
text(0.55, 0.88, 'Peak Ca', 'FontSize', 9, 'FontWeight', 'bold');
text(0.70, 0.88, 'Ca OK?', 'FontSize', 9, 'FontWeight', 'bold');
text(0.85, 0.88, 'Optimal', 'FontSize', 9, 'FontWeight', 'bold');

% Draw header line
line([0.05 0.95], [0.86 0.86], 'Color', 'k', 'LineWidth', 1);

% Table data - reduced spacing to fit all rows above the summary
row_spacing = 0.042;
start_y = 0.82;

for m = 1:n_Mg_hr
    y = start_y - (m-1)*row_spacing;
    
    % Highlight optimal rows with background
    if optimal(m)
        rectangle('Position', [0.05, y-0.015, 0.9, 0.035], ...
            'FaceColor', [0.9 0.95 0.9], 'EdgeColor', 'none');
    end
    
    text(0.10, y, sprintf('%.1f mM', Mg_highres(m)), 'FontSize', 8);
    text(0.25, y, sprintf('%d', hr_results.spikes(m)), 'FontSize', 8, 'HorizontalAlignment', 'center');
    text(0.40, y, sprintf('%.1f', hr_results.spike_red(m)), 'FontSize', 8, 'HorizontalAlignment', 'center');
    
    if ca_safe(m)
        ca_color = [0.2 0.6 0.2];
    else
        ca_color = [0.8 0.2 0.2];
    end
    text(0.55, y, sprintf('%.2f', hr_results.peak_Ca(m)), 'FontSize', 8, 'Color', ca_color, 'HorizontalAlignment', 'center');
    
    if ca_safe(m)
        text(0.70, y, 'YES', 'FontSize', 8, 'Color', [0.2 0.6 0.2], 'HorizontalAlignment', 'center');
    else
        text(0.70, y, 'NO', 'FontSize', 8, 'Color', [0.8 0.2 0.2], 'HorizontalAlignment', 'center');
    end
    
    if optimal(m)
        text(0.85, y, '✓', 'FontSize', 10, 'Color', [0.2 0.6 0.2], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

% Draw bottom line above summary
line([0.05 0.95], [0.12 0.12], 'Color', 'k', 'LineWidth', 1);

% Summary at bottom - well separated from data
if ~isempty(opt_idx)
    text(0.5, 0.07, sprintf('OPTIMAL RANGE: %.1f – %.1f mM', Mg_highres(opt_idx(1)), Mg_highres(opt_idx(end))), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2 0.5 0.2], 'HorizontalAlignment', 'center');
end

text(0.5, 0.02, sprintf('Criteria: Peak Ca < %.1f µM AND spike loss ≤ %d%%', Ca_toxic, max_spike_loss), ...
    'FontSize', 9, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');

sgtitle('Figure S2: Detailed Results', 'FontSize', 12, 'FontWeight', 'bold');

%% =========================================================================
%  SAVE FIGURES
%  =========================================================================
fprintf('Saving supplementary figures...\n');

save_figure(figS1, 'FigureS1_HighResolution');
save_figure(figS2, 'FigureS2_ResultsTable');

fprintf('Done!\n');

%% =========================================================================
%  PRINT SUMMARY
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  HIGH-RESOLUTION RESULTS (80 Hz)\n');
fprintf('=========================================================================\n\n');

fprintf('%-8s | %-8s | %-10s | %-10s | %-8s | %-8s\n', ...
    '[Mg²⁺]', 'Spikes', 'Spike Loss', 'Peak Ca', 'Ca OK?', 'Optimal');
fprintf('%s\n', repmat('-', 1, 60));

for m = 1:n_Mg_hr
    ca_str = 'NO'; opt_str = '';
    if ca_safe(m), ca_str = 'YES'; end
    if optimal(m), opt_str = '***'; end
    
    fprintf('%-8.1f | %-8d | %-10.1f | %-10.2f | %-8s | %-8s\n', ...
        Mg_highres(m), hr_results.spikes(m), hr_results.spike_red(m), ...
        hr_results.peak_Ca(m), ca_str, opt_str);
end

if ~isempty(opt_idx)
    fprintf('\n✅ OPTIMAL RANGE: %.1f – %.1f mM\n', Mg_highres(opt_idx(1)), Mg_highres(opt_idx(end)));
end

%% =========================================================================
%  FUNCTIONS
%  =========================================================================
function save_figure(fig, filename)
    % Save PNG
    try
        saveas(fig, [filename '.png']);
        fprintf('  Saved: %s.png\n', filename);
    catch ME
        warning('Could not save PNG: %s', ME.message);
    end
    
    % Save PDF with error handling
    try
        set(fig, 'PaperPositionMode', 'auto');
        fig_pos = fig.Position;
        set(fig, 'PaperUnits', 'points');
        set(fig, 'PaperSize', [fig_pos(3) fig_pos(4)]);
        print(fig, [filename '.pdf'], '-dpdf', '-vector');
        fprintf('  Saved: %s.pdf\n', filename);
    catch ME
        warning('Could not save PDF (file may be open or folder write-protected): %s', ME.message);
        fprintf('  PDF not saved - try closing the file if open, or save to a different folder.\n');
    end
end

function [V, Ca] = run_sim(t, dt, Glu, Mg_ext, p)
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
        
        Ca_influx = -p.k_Ca_NMDA * p.f_Ca * I_NMDA - p.k_Ca_CaL * I_CaL;
        Ca(i+1) = max(p.Ca_rest, Ca(i) + dt*(Ca_influx - (Ca(i)-p.Ca_rest)/p.tau_Ca));
        
        I_total = I_Na + I_Kdr + I_KA + I_CaL + I_KCa + I_L + I_AMPA + I_NMDA;
        V(i+1) = v + dt*(-I_total)/p.Cm;
    end
end
