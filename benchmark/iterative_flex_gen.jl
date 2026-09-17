using DataFrames, CSV, Plots, ExaModels, ExaModelsPower, MadNLP, MadNLPGPU, ExaPowerIO, CUDA

function find_shedding_limits(pg_solutions, filename, save_filename; max_samples=30, generator_order= [], best_order = false, save_folder = "saved_data/", shed_penalty=1e6)
    parsed_data = ExaPowerIO.parse_matpower(filename;)# library=:pglib)
    df = DataFrame(
            initial_solution = String[],
            nflex = Int[],
            obj = Float64[],
            total_pshed = Float64[],
            max_pshed = Float64[],
            total_qshed = Float64[],
            max_qshed = Float64[],
            flex_gen = []
    )

    gens = parsed_data.gen
    L_gen = length(parsed_data.gen)

    if isempty(generator_order)
        generator_order = reverse(sort(1:length(gens), by = i -> gens[i].pmax))
    end

    flex_generators = Int[]


    if L_gen ≤ max_samples
        gen_range = collect(0:L_gen)
    else
        gen_range = unique(round.(Int, range(0, L_gen, length=max_samples)))
    end

    gen_order_save = Int[]
    count_adds_save = Int[]
    added_values_save = Int[]
    
    for sol_tuple in pg_solutions


        
        if best_order
            generator_order, count_adds, added_values = find_best_ordering(filename, sol_tuple.pg_sol; c_shed = shed_penalty)
        end

        if sol_tuple.initial_method == "ACOPF tol 1e-4"
            gen_order_save = generator_order
            count_adds_save = count_adds
            added_values_save = added_values
        end


        for gen in gen_range
            println(sol_tuple.initial_method)
            println(save_filename)
            println(gen)
            j=0
            if gen == 1   
                flex_generators = [generator_order[1]]
            elseif gen == 0
                flex_generators = Int[]
            else
                flex_generators = generator_order[1:gen]
            end
            flex_gen = sum([parsed_data.gen[g].pmax for g in flex_generators]; init=0)
            while j < 3
                m_mod, v_mod, c_mod = build_augmented_opf(filename, flex_generators, shed_penalty, sol_tuple.pg_sol; backend = CUDABackend())
                result = madnlp(m_mod, tol=1e-8)
                p_sheds = Array(solution(result, v_mod.p_shed))
                q_sheds = Array(solution(result, v_mod.q_shed))
                #this needs to check if solver converged
                if result.status == MadNLP.SOLVE_SUCCEEDED
                    j = 4
                    push!(df, [sol_tuple.initial_method, length(flex_generators), result.objective, sum(abs.(p_sheds)), maximum(abs.(p_sheds)), sum(abs.(q_sheds)), maximum(abs.(q_sheds)), flex_gen])
                else 
                    j += 1
                end
            end
        end
    end

    CSV.write(save_folder*save_filename, df)
    return gen_order_save, count_adds_save, added_values_save
end

function test_ac_dc_combo(filename, case_name; max_samples=30, sort_order = :default, save_folder = "saved_data/", shed_penalty = 1e6)
    m_ac, v_ac, c_ac = ac_opf_model(filename; backend = CUDABackend())
    solved_limits = []

    result_tight = madnlp(m_ac, tol=1e-8)
    if result_tight.status == MadNLP.SOLVE_SUCCEEDED
        tight_tuple = (pg_sol=Array(solution(result_tight, v_ac.pg)), initial_method = "ACOPF tol 1e-8")
        push!(solved_limits, tight_tuple)
    else
        return nothing
    end
    
    result_loose = madnlp(m_ac, tol=1e-4)
    m_dc, v_dc, c_dc = dcopf_model(filename; backend=CUDABackend())
    result_dc = madnlp(m_dc, tol=1e-8)
    result_soc = solve_opf(filename, SOCWRPowerModel, Ipopt.Optimizer)
    result_qc = solve_opf(filename, QCRMPowerModel, Ipopt.Optimizer)
    loose_tuple = (pg_sol=Array(solution(result_loose, v_ac.pg)), initial_method = "ACOPF tol 1e-4")
    dc_tuple = (pg_sol=Array(solution(result_dc, v_dc.pg)), initial_method = "DCOPF tol 1e-8")
    socwr_tuple = (pg_sol=powermodels_pg_result(result_soc), initial_method = "SOCWR tol 1e-8")
    qcrm_tuple = (pg_sol=powermodels_pg_result(result_qc), initial_method = "QCRM tol 1e-8")
    #sdpwrm_tuple = (pg_sol=powermodels_pg_result(filename, SDPWRMPowerModel), initial_method = "SDPWRM tol 1e-8")


    save_filename = save_folder*case_name

    best_order=false

    if sort_order == :ac
        ordered_gens = find_ac_gen_order(tight_tuple.pg_sol)
        case_name = case_name*"_ac_ordered"
    elseif sort_order == :best
        best_order = true
        ordered_gens = []
        case_name = case_name*"_best_ordered"
    else
        ordered_gens = []
    end

    
    if result_loose.status == MadNLP.SOLVE_SUCCEEDED
        push!(solved_limits, loose_tuple)
    end
    if result_dc.status == MadNLP.SOLVE_SUCCEEDED
        push!(solved_limits, dc_tuple)
    end
    if result_soc["termination_status"] == LOCALLY_SOLVED
        push!(solved_limits, socwr_tuple)
    end
    if result_qc["termination_status"] == LOCALLY_SOLVED
        push!(solved_limits, qcrm_tuple)
    end

    gen_order_save, count_adds_save, added_values = find_shedding_limits(solved_limits, filename, case_name*"_samples"*string(max_samples); max_samples=max_samples, generator_order=ordered_gens, best_order=best_order, save_folder=save_folder, shed_penalty=shed_penalty)
    group_plot(case_name, max_samples, save_folder)
    return gen_order_save, count_adds_save, added_values
end

function group_plot(case_name, max_samples, save_folder)
    saved_filename = save_folder*case_name*"_samples"*string(max_samples)
    plot_nflex_vs_y(saved_filename, "obj", save_folder*case_name*"_samples"*string(max_samples)*"_obj"; title="", ylabel = "Objective")
    plot_nflex_vs_y(saved_filename, "obj", save_folder*case_name*"_samples"*string(max_samples)*"_obj"*"_short"; title="", tight_bounds=true, ylabel = "Objective")
    plot_nflex_vs_y(saved_filename, "total_pshed", save_folder*case_name*"_samples"*string(max_samples)*"_p_tot"; title="", ylabel = "Total Active Power Load Shed")
    plot_nflex_vs_y(saved_filename, "max_pshed", save_folder*case_name*"_samples"*string(max_samples)*"_p_max"; title="", ylabel = "Maximum Active Power Load Shed")
end

function plot_nflex_vs_y(
    csvfile::AbstractString,
    ycol::String,
    savepath::AbstractString;
    xlabel::AbstractString = "Number of Adjustible Generators",
    ylabel::AbstractString = ycol,
    title::AbstractString = "",
    tight_bounds = false
)
    gr()              # reinitialize backend
    Plots.closeall() 

    df = CSV.read(csvfile, DataFrame)


    @assert "initial_solution" in names(df)
    @assert "nflex" in names(df)
    @assert ycol in names(df)



    p = plot(
        xlabel = xlabel,
        ylabel = ylabel,
        yscale = :log10,
        title = title,
        legend = :best,
    )



    linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    markers    = [:circle, :square, :diamond, :utriangle, :dtriangle, :cross]
    alphas     = range(0.9, 0.5; length=length(unique(df.initial_solution)))

    #if tight_bounds
    #    baseline_row = findfirst(df[!, "initial_solution"] .== "ACOPF tol 1e-8")
    #    if baseline_row !== nothing
    #        baseline = df[baseline_row, ycol]
    #        ylims = (0.99 * baseline, 1.05 * baseline)
    #        plot!(p, ylim = ylims, yscale = :linear,)
    #    end
    #end


    for (i, sol) in enumerate(unique(df.initial_solution))
        sub = df[df[!, "initial_solution"] .== sol, :]
        if tight_bounds
            sub = sub[sub[!, "max_pshed"] .< 1e-7, :]
            plot!(p, yscale = :linear,)
        end
        sort!(sub, "nflex")

        plot!(
            p,
            sub[!, "nflex"],
            sub[!, ycol],
            label     = sol,
            linestyle = linestyles[mod1(i, length(linestyles))],
            marker    = markers[mod1(i, length(markers))],
            alpha     = alphas[i],
            linewidth = 2,
            markersize = 4,
        )
    end


    savefig(p, savepath*".png")
end

function find_ac_gen_order(pg_sol)
    pg_sol = Array(pg_sol)
    return sortperm(pg_sol, rev = true)
end