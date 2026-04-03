function outFile = export_case118_ddvpp_to_excel(outFile)
% EXPORT_CASE118_DDVPP_TO_EXCEL
% Export MATPOWER case118 + DDVPP setup parameters to a single Excel file.
%
% Usage:
%   export_case118_ddvpp_to_excel
%   export_case118_ddvpp_to_excel('case118_ddvpp_parameters.xlsx')
%
% Requirements:
%   1) MATPOWER is on the MATLAB path, including case118 and makeBdc.
%   2) setup_ddvpp_case118_ddvpp.m is on the MATLAB path.
%   3) Optional: ddvpp_build_model_data.m on the MATLAB path to export
%      reduced-model matrices J, L, F, and Dload_equiv.
%
% Output:
%   An Excel workbook containing raw MATPOWER data, DDVPP settings, and
%   optional reduced-model matrices for later direct reuse.

    if nargin < 1 || isempty(outFile)
        outFile = 'case118_ddvpp_parameters.xlsx';
    end

    assert(exist('setup_ddvpp_case118_ddvpp', 'file') == 2, ...
        'setup_ddvpp_case118_ddvpp.m not found on MATLAB path.');
    assert(exist('case118', 'file') == 2, ...
        'case118.m not found on MATLAB path. Please add MATPOWER to the path.');

    mpc = setup_ddvpp_case118_ddvpp();

    if exist(outFile, 'file') == 2
        delete(outFile);
    end

    %% Summary sheet
    gen_buses = mpc.gen(:,1);
    all_buses = mpc.bus(:,1);
    load_buses = setdiff(all_buses, gen_buses, 'stable');
    if isfield(mpc, 'ddvpp_bounds')
        ibr_buses = mpc.ddvpp_bounds(:,1);
    else
        ibr_buses = [];
    end

    summary = {
        'Field', 'Value';
        'baseMVA', mpc.baseMVA;
        'num_bus_rows', size(mpc.bus,1);
        'num_gen_rows', size(mpc.gen,1);
        'num_branch_rows', size(mpc.branch,1);
        'num_load_side_buses', numel(load_buses);
        'num_IBR_buses', numel(ibr_buses);
        'has_gencost', isfield(mpc,'gencost');
        'has_dyn', isfield(mpc,'dyn');
        'has_ddvpp_bounds', isfield(mpc,'ddvpp_bounds');
        'has_mu_load', isfield(mpc,'mu_load');
        'has_security', isfield(mpc,'security');
        'has_admm', isfield(mpc,'admm')
        };
    writecell(summary, outFile, 'Sheet', 'summary', 'Range', 'A1');

    %% Raw MATPOWER sheets
    write_numeric_sheet(outFile, 'bus', mpc.bus, bus_headers(size(mpc.bus,2)));
    write_numeric_sheet(outFile, 'gen', mpc.gen, gen_headers(size(mpc.gen,2)));
    write_numeric_sheet(outFile, 'branch', mpc.branch, branch_headers(size(mpc.branch,2)));

    if isfield(mpc, 'gencost')
        write_numeric_sheet(outFile, 'gencost', mpc.gencost, generic_headers('gencost_col_', size(mpc.gencost,2)));
    end

    %% DDVPP setup sheets
    if isfield(mpc, 'dyn')
        dyn_hdr = {'bus_id','type_1SG_2IBR','m0','d0','k','tau','gamma'};
        write_numeric_sheet(outFile, 'dyn', mpc.dyn, dyn_hdr);
    end

    if isfield(mpc, 'ddvpp_bounds')
        bnd_hdr = {'bus_id','m_max','d_max','c_m','c_d'};
        write_numeric_sheet(outFile, 'ddvpp_bounds', mpc.ddvpp_bounds, bnd_hdr);
    end

    if isfield(mpc, 'mu_load')
        mu_tbl = [load_buses(:), mpc.mu_load(:)];
        write_numeric_sheet(outFile, 'mu_load', mu_tbl, {'load_bus_id','mu_load'});
    end

    if isfield(mpc, 'security')
        write_struct_sheet(outFile, 'security', mpc.security);
    end

    if isfield(mpc, 'admm')
        write_struct_sheet(outFile, 'admm', mpc.admm);
    end

    %% Bus classification sheet
    bus_sets = zeros(numel(all_buses), 4);
    bus_sets(:,1) = all_buses(:);
    bus_sets(:,2) = ismember(all_buses, gen_buses);
    bus_sets(:,3) = ismember(all_buses, load_buses);
    bus_sets(:,4) = ismember(all_buses, ibr_buses);
    write_numeric_sheet(outFile, 'bus_sets', bus_sets, {'bus_id','is_gen_bus','is_load_bus','is_ibr_bus'});

    %% Optional reduced-model export
    if exist('ddvpp_build_model_data', 'file') == 2 && exist('makeBdc', 'file') == 2
        try
            data = ddvpp_build_model_data(mpc);
            meta = {
                'Field', 'Value';
                'ngen', data.ngen;
                'nload', data.nload;
                'nbus', data.nbus;
                'f_base', data.f_base
                };
            writecell(meta, outFile, 'Sheet', 'reduced_meta', 'Range', 'A1');

            write_labeled_matrix_sheet(outFile, 'J', data.J, data.gen_buses, data.gen_buses, 'row_bus_id', 'col_bus_');
            write_labeled_matrix_sheet(outFile, 'L', data.L, data.gen_buses, data.load_buses, 'row_gen_bus_id', 'col_load_bus_');
            write_labeled_matrix_sheet(outFile, 'F', data.F, data.load_buses, data.gen_buses, 'row_load_bus_id', 'col_gen_bus_');
            write_labeled_matrix_sheet(outFile, 'Dload_equiv', data.Dload_equiv, data.gen_buses, data.gen_buses, 'row_bus_id', 'col_bus_');
        catch ME
            warn = {
                'Warning';
                ['Reduced-model export skipped: ' ME.message]
                };
            writecell(warn, outFile, 'Sheet', 'reduced_meta', 'Range', 'A8');
        end
    else
        note = {
            'Warning';
            'Reduced-model export skipped because ddvpp_build_model_data.m or makeBdc.m is missing from MATLAB path.'
            };
        writecell(note, outFile, 'Sheet', 'reduced_meta', 'Range', 'A1');
    end

    fprintf('Excel export completed: %s\n', outFile);
end

function write_numeric_sheet(outFile, sheetName, M, headers)
    if isempty(M)
        writecell(headers, outFile, 'Sheet', sheetName, 'Range', 'A1');
        return;
    end
    T = array2table(M, 'VariableNames', matlab.lang.makeValidName(headers, 'ReplacementStyle', 'delete'));
    writetable(T, outFile, 'Sheet', sheetName, 'WriteMode', 'overwrite');
end

function write_struct_sheet(outFile, sheetName, s)
    f = fieldnames(s);
    vals = cell(numel(f), 1);
    for i = 1:numel(f)
        v = s.(f{i});
        if islogical(v)
            vals{i} = double(v);
        elseif isnumeric(v) && isscalar(v)
            vals{i} = v;
        elseif ischar(v) || isstring(v)
            vals{i} = char(string(v));
        else
            vals{i} = evalc('disp(v)');
            vals{i} = strtrim(vals{i});
        end
    end
    C = [{'field','value'}; [f, vals]];
    writecell(C, outFile, 'Sheet', sheetName, 'Range', 'A1');
end

function write_labeled_matrix_sheet(outFile, sheetName, M, rowIDs, colIDs, firstColHeader, colPrefix)
    colHeaders = cell(1, numel(colIDs) + 1);
    colHeaders{1} = firstColHeader;
    for j = 1:numel(colIDs)
        colHeaders{j+1} = sprintf('%s%d', colPrefix, colIDs(j));
    end
    C = cell(size(M,1) + 1, size(M,2) + 1);
    C(1,:) = colHeaders;
    for i = 1:size(M,1)
        C{i+1,1} = rowIDs(i);
        for j = 1:size(M,2)
            C{i+1,j+1} = M(i,j);
        end
    end
    writecell(C, outFile, 'Sheet', sheetName, 'Range', 'A1');
end

function hdr = bus_headers(n)
    base = {'BUS_I','BUS_TYPE','PD','QD','GS','BS','BUS_AREA','VM','VA','BASE_KV','ZONE','VMAX','VMIN'};
    hdr = fill_headers(base, n, 'bus_extra_');
end

function hdr = gen_headers(n)
    base = {'GEN_BUS','PG','QG','QMAX','QMIN','VG','MBASE','GEN_STATUS','PMAX','PMIN', ...
            'PC1','PC2','QC1MIN','QC1MAX','QC2MIN','QC2MAX','RAMP_AGC','RAMP_10','RAMP_30','RAMP_Q','APF'};
    hdr = fill_headers(base, n, 'gen_extra_');
end

function hdr = branch_headers(n)
    base = {'F_BUS','T_BUS','BR_R','BR_X','BR_B','RATE_A','RATE_B','RATE_C','TAP','SHIFT','BR_STATUS','ANGMIN','ANGMAX'};
    hdr = fill_headers(base, n, 'branch_extra_');
end

function hdr = generic_headers(prefix, n)
    hdr = cell(1,n);
    for i = 1:n
        hdr{i} = sprintf('%s%d', prefix, i);
    end
end

function hdr = fill_headers(base, n, prefix)
    hdr = cell(1,n);
    nBase = min(numel(base), n);
    hdr(1:nBase) = base(1:nBase);
    for i = nBase+1:n
        hdr{i} = sprintf('%s%d', prefix, i - nBase);
    end
end
