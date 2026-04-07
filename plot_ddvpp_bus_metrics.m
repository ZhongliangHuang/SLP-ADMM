function figs = plot_ddvpp_bus_metrics(eval_out, save_dir, topk_curve)
%PLOT_DDVPP_BUS_METRICS
% Plot bus-wise RoCoF, nadir, nadir time, and frequency-response curves
% from evaluate_ddvpp_frequency_response output.
%
% Inputs
%   eval_out    : output struct from evaluate_ddvpp_frequency_response(...)
%   save_dir    : optional folder path for saving png figures
%   topk_curve  : optional number of worst buses to show in detailed curves
%
% Outputs
%   figs        : struct containing figure handles
%
% The function is written to be robust to different versions of the
% evaluator output structure.

    if nargin < 2 || isempty(save_dir)
        save_dir = '';
    end
    if nargin < 3 || isempty(topk_curve)
        topk_curve = 10;
    end

    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
    end

    %% --------------------------------------------------------------------
    %% Extract per-bus metrics table
    %% --------------------------------------------------------------------
    T = local_get_metrics_table(eval_out);

    bus = T.bus;
    rocof = T.rocof0_hz_per_s;
    nadir = T.nadir_hz;
    nadir_time = T.nadir_time_s;
    abs_rocof = abs(rocof);

    %% --------------------------------------------------------------------
    %% Extract time axis and load-side frequency trajectories in Hz
    %% --------------------------------------------------------------------
    t = local_get_time_vector(eval_out);
    resp_load_hz = local_get_resp_load_hz(eval_out);

    %% --------------------------------------------------------------------
    %% Identify worst buses
    %% --------------------------------------------------------------------
    [~, idx_worst_rocof] = max(abs_rocof);
    [~, idx_worst_nadir] = min(nadir);
    [~, idx_latest_nadir] = max(nadir_time);

    [~, order_rocof] = sort(abs_rocof, 'descend');
    [~, order_nadir] = sort(nadir, 'ascend');

    topk_curve = min(topk_curve, numel(bus));
    topk_bar = min(10, numel(bus));

    idx_top_rocof_curve = order_rocof(1:topk_curve);
    idx_top_nadir_curve = order_nadir(1:topk_curve);

    idx_top_rocof_bar = order_rocof(1:topk_bar);
    idx_top_nadir_bar = order_nadir(1:topk_bar);

    figs = struct();

    %% --------------------------------------------------------------------
    %% Figure 1: Bus-wise RoCoF
    %% --------------------------------------------------------------------
    figs.rocof = figure('Name', 'Bus-wise RoCoF', 'Color', 'w');
    plot(bus, rocof, 'o-', 'LineWidth', 1.0, 'MarkerSize', 4);
    hold on;
    plot(bus(idx_worst_rocof), rocof(idx_worst_rocof), 's', ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    xlabel('Bus index');
    ylabel('RoCoF at t = 0^+ (Hz/s)');
    title('Bus-wise initial RoCoF');
    text(bus(idx_worst_rocof), rocof(idx_worst_rocof), ...
        sprintf('  worst bus = %d, %.4f Hz/s', bus(idx_worst_rocof), rocof(idx_worst_rocof)), ...
        'FontSize', 10, 'VerticalAlignment', 'bottom');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.rocof, fullfile(save_dir, 'bus_rocof.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 2: Bus-wise nadir
    %% --------------------------------------------------------------------
    figs.nadir = figure('Name', 'Bus-wise nadir', 'Color', 'w');
    plot(bus, nadir, 'o-', 'LineWidth', 1.0, 'MarkerSize', 4);
    hold on;
    plot(bus(idx_worst_nadir), nadir(idx_worst_nadir), 's', ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    xlabel('Bus index');
    ylabel('Frequency nadir (Hz)');
    title('Bus-wise frequency nadir');
    text(bus(idx_worst_nadir), nadir(idx_worst_nadir), ...
        sprintf('  worst bus = %d, %.4f Hz', bus(idx_worst_nadir), nadir(idx_worst_nadir)), ...
        'FontSize', 10, 'VerticalAlignment', 'top');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.nadir, fullfile(save_dir, 'bus_nadir.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 3: Bus-wise nadir time
    %% --------------------------------------------------------------------
    figs.nadir_time = figure('Name', 'Bus-wise nadir time', 'Color', 'w');
    plot(bus, nadir_time, 'o-', 'LineWidth', 1.0, 'MarkerSize', 4);
    hold on;
    plot(bus(idx_latest_nadir), nadir_time(idx_latest_nadir), 's', ...
        'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    xlabel('Bus index');
    ylabel('Time to nadir (s)');
    title('Bus-wise nadir arrival time');
    text(bus(idx_latest_nadir), nadir_time(idx_latest_nadir), ...
        sprintf('  latest nadir bus = %d, %.4f s', bus(idx_latest_nadir), nadir_time(idx_latest_nadir)), ...
        'FontSize', 10, 'VerticalAlignment', 'bottom');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.nadir_time, fullfile(save_dir, 'bus_nadir_time.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 4: All load-bus frequency-response curves
    %% --------------------------------------------------------------------
    figs.all_curves = figure('Name', 'All bus frequency-response curves', 'Color', 'w');
    plot(t, resp_load_hz.', 'LineWidth', 0.8);
    grid on;
    xlabel('Time (s)');
    ylabel('Frequency deviation (Hz)');
    title('All load-side bus frequency responses');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.all_curves, fullfile(save_dir, 'all_bus_frequency_curves.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 5: Top-k worst nadir bus curves
    %% --------------------------------------------------------------------
    figs.top_nadir_curves = figure('Name', 'Top worst nadir bus curves', 'Color', 'w');
    plot(t, resp_load_hz(idx_top_nadir_curve,:).', 'LineWidth', 1.2);
    grid on;
    xlabel('Time (s)');
    ylabel('Frequency deviation (Hz)');
    title(sprintf('Top-%d worst buses by nadir', topk_curve));

    legend_labels = arrayfun(@(x) sprintf('Bus %d', x), bus(idx_top_nadir_curve), 'UniformOutput', false);
    legend(legend_labels, 'Location', 'bestoutside');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.top_nadir_curves, fullfile(save_dir, 'top_nadir_bus_curves.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 6: Top-k worst |RoCoF| bus curves
    %% --------------------------------------------------------------------
    figs.top_rocof_curves = figure('Name', 'Top worst RoCoF bus curves', 'Color', 'w');
    plot(t, resp_load_hz(idx_top_rocof_curve,:).', 'LineWidth', 1.2);
    grid on;
    xlabel('Time (s)');
    ylabel('Frequency deviation (Hz)');
    title(sprintf('Top-%d worst buses by |RoCoF|', topk_curve));

    legend_labels = arrayfun(@(x) sprintf('Bus %d', x), bus(idx_top_rocof_curve), 'UniformOutput', false);
    legend(legend_labels, 'Location', 'bestoutside');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.top_rocof_curves, fullfile(save_dir, 'top_rocof_bus_curves.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 7: Top-10 worst |RoCoF|
    %% --------------------------------------------------------------------
    figs.top_rocof = figure('Name', 'Top-10 worst RoCoF buses', 'Color', 'w');
    bar(categorical(string(bus(idx_top_rocof_bar))), rocof(idx_top_rocof_bar));
    grid on;
    xlabel('Bus index');
    ylabel('RoCoF at t = 0^+ (Hz/s)');
    title('Top worst buses by |RoCoF|');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.top_rocof, fullfile(save_dir, 'top10_rocof.png'), 'Resolution', 300);
    % end

    %% --------------------------------------------------------------------
    %% Figure 8: Top-10 worst nadir
    %% --------------------------------------------------------------------
    figs.top_nadir = figure('Name', 'Top-10 worst nadir buses', 'Color', 'w');
    bar(categorical(string(bus(idx_top_nadir_bar))), nadir(idx_top_nadir_bar));
    grid on;
    xlabel('Bus index');
    ylabel('Frequency nadir (Hz)');
    title('Top worst buses by nadir');

    % if ~isempty(save_dir)
    %     exportgraphics(figs.top_nadir, fullfile(save_dir, 'top10_nadir.png'), 'Resolution', 300);
    % end
end

%% ========================================================================
function T = local_get_metrics_table(eval_out)
% Robustly extract the per-bus metrics table

    if isfield(eval_out, 'metrics')
        if istable(eval_out.metrics)
            T = eval_out.metrics;
            return;
        elseif isstruct(eval_out.metrics)
            if isfield(eval_out.metrics, 'per_bus') && istable(eval_out.metrics.per_bus)
                T = eval_out.metrics.per_bus;
                return;
            end
        end
    end

    error(['Cannot find per-bus metrics table. Expected one of: ', ...
           'eval_out.metrics.per_bus or eval_out.metrics']);
end

%% ========================================================================
function t = local_get_time_vector(eval_out)
% Robustly extract time vector

    if isfield(eval_out, 'response') && isfield(eval_out.response, 't')
        t = eval_out.response.t;
        return;
    end

    if isfield(eval_out, 'result') && isfield(eval_out.result, 't')
        t = eval_out.result.t;
        return;
    end

    error('Cannot find time vector t in eval_out.response or eval_out.result.');
end

%% ========================================================================
function resp_load_hz = local_get_resp_load_hz(eval_out)
% Robustly extract load-side frequency-response trajectories in Hz

    % Case 1: already stored in Hz
    if isfield(eval_out, 'response') && isfield(eval_out.response, 'resp_load_hz')
        resp_load_hz = eval_out.response.resp_load_hz;
        return;
    end

    if isfield(eval_out, 'result') && isfield(eval_out.result, 'resp_load_hz')
        resp_load_hz = eval_out.result.resp_load_hz;
        return;
    end

    % Case 2: stored in p.u., convert using base frequency
    if isfield(eval_out, 'response') && isfield(eval_out.response, 'resp_load_pu')
        resp_load_pu = eval_out.response.resp_load_pu;
    elseif isfield(eval_out, 'result') && isfield(eval_out.result, 'resp_load_pu')
        resp_load_pu = eval_out.result.resp_load_pu;
    else
        error(['Cannot find load-side frequency trajectories. Expected one of: ', ...
               'resp_load_hz or resp_load_pu']);
    end

    f0 = local_get_base_frequency(eval_out);
    resp_load_hz = resp_load_pu * f0;
end

%% ========================================================================
function f0 = local_get_base_frequency(eval_out)
% Robustly extract base frequency from several possible locations

    if isfield(eval_out, 'cfg') && isfield(eval_out.cfg, 'base_frequency_hz')
        f0 = eval_out.cfg.base_frequency_hz;
        return;
    end

    if isfield(eval_out, 'model') && isfield(eval_out.model, 'cfg') ...
            && isfield(eval_out.model.cfg, 'base_frequency_hz')
        f0 = eval_out.model.cfg.base_frequency_hz;
        return;
    end

    if isfield(eval_out, 'model') && isfield(eval_out.model, 'mpc') ...
            && isfield(eval_out.model.mpc, 'userdata') ...
            && isfield(eval_out.model.mpc.userdata, 'dynamic') ...
            && isfield(eval_out.model.mpc.userdata.dynamic, 'base_frequency_hz')
        f0 = eval_out.model.mpc.userdata.dynamic.base_frequency_hz;
        return;
    end

    if isfield(eval_out, 'mpc') && isfield(eval_out.mpc, 'userdata') ...
            && isfield(eval_out.mpc.userdata, 'dynamic') ...
            && isfield(eval_out.mpc.userdata.dynamic, 'base_frequency_hz')
        f0 = eval_out.mpc.userdata.dynamic.base_frequency_hz;
        return;
    end

    error(['base_frequency_hz not found. Expected one of: ', ...
           'eval_out.cfg.base_frequency_hz, ', ...
           'eval_out.model.cfg.base_frequency_hz, ', ...
           'eval_out.model.mpc.userdata.dynamic.base_frequency_hz, ', ...
           'or eval_out.mpc.userdata.dynamic.base_frequency_hz']);
end