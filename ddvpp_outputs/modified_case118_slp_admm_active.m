function mpc = modified_case118_slp_admm_active()
%MODIFIED_CASE118_SLP_ADMM_ACTIVE
% -------------------------------------------------------------------------
% 这是一个“让 SLP + ADMM 真正活跃起来”的 case 包装器。
%
% 为什么需要它：
% 当前仓库里的 modified_case118.m 默认安全限值较松：
%   - nadir_limit_hz = 0.8
%   - qss_limit_hz   = 0.2
% 而你前面的基线结果大约是：
%   - worst nadir ≈ -0.174 Hz
%   - worst QSS   ≈ 0.090 Hz
% 因此在这个 case 上，nadir / QSS 约束天然不激活，
% SLP 外层会长期得到 req Delta sigma = 0、req sum(Delta D) = 0，
% 最终只剩 RoCoF 在推着算法走。
%
% 这个包装器做三件事：
%   1) 收紧 nadir / QSS 限值，使 sigma-left-shift 与 Delta D 进入活跃集；
%   2) 把 row-1 扰动改到更脆弱的北部弱耦合区域附近，让空间差异更明显；
%   3) 给关键 IBR 节点更大的 M / D headroom 和更低的局部成本，
%      让 ADMM 的“有差异分配”更容易体现出来。
%
% 推荐用法：
%   p = struct();
%   p.case_function = 'modified_case118_slp_admm_active';
%   p.disturbance_id = 1;
%   result = run_ddvpp_slp_admm_design_v2(p);
%
% -------------------------------------------------------------------------

mpc = modified_case118();
gtab = mpc.userdata.ddvpp.gen_dynamic_table;

% -------------------------------------------------------------------------
% 1) 收紧安全约束，让 nadir / QSS 真正进入活跃集
% -------------------------------------------------------------------------
mpc.userdata.ddvpp.security_limits.nadir_limit_hz = 0.22;
mpc.userdata.ddvpp.security_limits.qss_limit_hz = 0.12;
mpc.userdata.ddvpp.security_limits.rocof_limit_hz_per_s = 1.0;

% 增加外层/内层次数上限，给更严格 case 足够迭代空间
mpc.userdata.ddvpp.security_limits.max_iterations_slp = 25;
mpc.userdata.ddvpp.security_limits.max_iterations_admm = 250;

% -------------------------------------------------------------------------
% 2) 调整 row-1 扰动：把默认测试场景移到更容易触发空间差异的位置
% -------------------------------------------------------------------------
% 当前原始 case 的 disturbance_set 为：
%   row1 bus=50, 1200 MW
%   row2 bus=49, 1000 MW
%   row3 bus=37,  800 MW
%   row4 bus=8,   800 MW
% 这里将 row1 改为 bus 49 的更强负荷扰动，优先观察 RoCoF 与局部 nadir。

baseMVA = mpc.baseMVA;
dset = mpc.userdata.ddvpp.disturbance_set;
dset.bus(1) = 49;
dset.deltaP_mw(1) = 1400.0;
dset.deltaP_pu(1) = dset.deltaP_mw(1) / baseMVA;
dset.weight(1) = 1.00;
mpc.userdata.ddvpp.disturbance_set = dset;

% 同步 evaluator 默认展示字段，便于不传参数时直接读 case 默认值
mpc.userdata.dynamic.disturbance_load_bus = 49;
mpc.userdata.dynamic.disturbance_mw = 1400.0;

% -------------------------------------------------------------------------
% 3) 给关键 IBR 节点更多 headroom，并且降低其局部成本
% -------------------------------------------------------------------------
% 这些节点本来就在当前 case 中被建模成 IBR 热点：
% [46 49 54 55 56] 处在你前面日志最敏感的区域附近。

critical_buses = [46 49 54 55 56];
mask_crit = ismember(gtab.host_bus, critical_buses) & gtab.is_controllable;
mask_other = gtab.is_controllable & ~mask_crit;

% 对关键节点放宽惯量/阻尼上界，避免只靠统一小步长慢慢推。
gtab.m_max(mask_crit) = 4.0 .* gtab.M0(mask_crit);
gtab.d_max(mask_crit) = 4.0 .* gtab.D0(mask_crit);

% 对非关键节点保持适度 headroom，防止“大家都一样加 M”。
gtab.m_max(mask_other) = 2.0 .* gtab.M0(mask_other);
gtab.d_max(mask_other) = 2.5 .* gtab.D0(mask_other);

% 给关键节点更低局部成本，让优化自然向高杠杆位置聚焦。
gtab.local_cost_quad(mask_crit) = 0.5 .* gtab.local_cost_quad(mask_crit);
gtab.local_cost_lin(mask_crit)  = 0.5 .* gtab.local_cost_lin(mask_crit);

% 对非关键节点略提高成本，增强差异化分配。
gtab.local_cost_quad(mask_other) = 1.2 .* gtab.local_cost_quad(mask_other);
gtab.local_cost_lin(mask_other)  = 1.2 .* gtab.local_cost_lin(mask_other);

mpc.userdata.ddvpp.gen_dynamic_table = gtab;

% 保留说明，方便以后回看 case 是怎么设计的。
notes = mpc.userdata.ddvpp.optimization_notes;
notes{end+1} = 'active wrapper: tightened nadir/qss limits to activate sigma and Delta-D constraints.';
notes{end+1} = 'active wrapper: row-1 disturbance moved to bus 49 and increased to 1400 MW.';
notes{end+1} = 'active wrapper: critical IBR buses [46 49 54 55 56] get larger M/D headroom and lower cost.';
mpc.userdata.ddvpp.optimization_notes = notes;
end
