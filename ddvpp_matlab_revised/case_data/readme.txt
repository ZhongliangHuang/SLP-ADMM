文件说明
1. export_case118_ddvpp_to_excel.m
   用于把 MATPOWER 的 case118.m 与 setup_ddvpp_case118_ddvpp.m 生成的系统参数导出到一个 Excel 工作簿中。
2. 本说明文件

运行前准备
1. 把 MATPOWER 加到 MATLAB 路径中，确保 case118.m 和 makeBdc.m 可用。
2. 把 setup_ddvpp_case118_ddvpp.m 放到 MATLAB 路径中。
3. 如果你还想一起导出降阶模型矩阵 J、L、F、Dload_equiv，请把 ddvpp_build_model_data.m 也放到 MATLAB 路径中。

推荐放在同一文件夹下的文件
- export_case118_ddvpp_to_excel.m
- setup_ddvpp_case118_ddvpp.m
- ddvpp_build_model_data.m
- 你的 MATPOWER 路径

使用方法
在 MATLAB 命令行运行：

export_case118_ddvpp_to_excel

或者指定输出文件名：

export_case118_ddvpp_to_excel('case118_ddvpp_parameters.xlsx')

生成的 Excel 工作簿包含这些工作表
1. summary
   基本统计信息。
2. bus
   case118 的 bus 矩阵。
3. gen
   case118 的 gen 矩阵。
4. branch
   case118 的 branch 矩阵。
5. gencost
   如果 mpc.gencost 存在，则一并导出。
6. dyn
   setup_ddvpp_case118_ddvpp.m 中的动态参数表。
7. ddvpp_bounds
   DDVPP 优化边界和成本参数。
8. mu_load
   负荷阻尼分布。
9. security
   安全约束参数。
10. admm
    ADMM 参数。
11. bus_sets
    每个母线是否为发电母线、负荷侧母线、IBR 母线。
12. reduced_meta
    降阶模型的维度信息和告警信息。
13. J
14. L
15. F
16. Dload_equiv
    如果 ddvpp_build_model_data.m 和 makeBdc.m 可用，则这些降阶矩阵会一并导出。

说明
1. 这个脚本不会修改原始 case118 数据，只会读取 setup_ddvpp_case118_ddvpp.m 返回的 mpc 结构并导出。
2. 如果 reduced_meta 中出现警告，一般表示 MATLAB 路径里缺少 ddvpp_build_model_data.m 或 makeBdc.m。
3. 如果你后续希望从 Excel 再读回 MATLAB，建议优先保留 bus、gen、branch、dyn、ddvpp_bounds、mu_load、security、admm 这些表。

建议后续扩展
1. 你可以再写一个反向读取脚本，例如 import_case118_ddvpp_from_excel.m。
2. 如果你想保留每次试验的参数快照，可以把输出文件名改成带日期的版本，例如：
   export_case118_ddvpp_to_excel(['case118_ddvpp_' datestr(now,'yyyymmdd_HHMMSS') '.xlsx'])
