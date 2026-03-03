using DataFrames, CSV, Plots, ExaModels, ExaModelsPower, MadNLP, MadNLPGPU, ExaPowerIO, CUDA

function find_shedding_limits(pg_solutions, filename, save_filename; max_samples=30, generator_order= [], best_order = false)
    parsed_data = ExaPowerIO.parse_matpower(filename;)# library=:pglib)
    df = DataFrame(
            initial_solution = String[],
            nflex = Int[],
            obj = Float64[],
            total_pshed = Float64[],
            max_pshed = Float64[],
            total_qshed = Float64[],
            max_qshed = Float64[]
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

    
    for sol_tuple in pg_solutions
        if sol_tuple.initial_method != "ACOPF tol 1e-4"
            continue
        end
        println(sol_tuple)
        if best_order
            generator_order = find_best_ordering(filename, sol_tuple.pg_sol)
        end


        skip = "
        for gen in gen_range
            j=0
            if gen == 1   
                flex_generators = [generator_order[1]]
            elseif gen == 0
                flex_generators = Int[]
            else
                flex_generators = generator_order[1:gen]
            end
            while j < 3
            
                println(length(flex_generators))
                println(sol_tuple.initial_method)
                println(gen)
                m_mod, v_mod, c_mod = build_augmented_opf(filename, flex_generators, 1e6, sol_tuple.pg_sol; backend = CUDABackend())
                result = madnlp(m_mod, tol=1e-8)
                p_sheds = Array(solution(result, v_mod.p_shed))
                q_sheds = Array(solution(result, v_mod.q_shed))
                #this needs to check if solver converged
                if result.status == MadNLP.SOLVE_SUCCEEDED
                    j = 4
                    push!(df, [sol_tuple.initial_method, length(flex_generators), result.objective, sum(abs.(p_sheds)), maximum(abs.(p_sheds)), sum(abs.(q_sheds)), maximum(abs.(q_sheds))])
                else 
                    j += 1
                end
            end
        end"
    end

    #CSV.write("data/"*save_filename, df)
end

function test_ac_dc_combo(filename, case_name; max_samples=30, sort_order = :default)
    m_ac, v_ac, c_ac = ac_opf_model(filename; backend = CUDABackend())
    result_loose = madnlp(m_ac, tol=1e-4)
    result_tight = madnlp(m_ac, tol=1e-8)
    m_dc, v_dc, c_dc = dcopf_model(filename; backend=CUDABackend())
    result_dc = madnlp(m_dc, tol=1e-8)
    loose_tuple = (pg_sol=Array(solution(result_loose, v_ac.pg)), initial_method = "ACOPF tol 1e-4")
    tight_tuple = (pg_sol=Array(solution(result_tight, v_ac.pg)), initial_method = "ACOPF tol 1e-8")
    dc_tuple = (pg_sol=Array(solution(result_dc, v_dc.pg)), initial_method = "DCOPF tol 1e-8")
    socwr_tuple = (pg_sol=powermodels_pg_result(filename, SOCWRPowerModel), initial_method = "SOCWR tol 1e-8")
    qcrm_tuple = (pg_sol=powermodels_pg_result(filename, QCRMPowerModel), initial_method = "QCRM tol 1e-8")
    #sdpwrm_tuple = (pg_sol=powermodels_pg_result(filename, SDPWRMPowerModel), initial_method = "SDPWRM tol 1e-8")


    save_filename = "data/"*case_name

    best_order=false

    case_name = case_name*"_samples"*string(max_samples)

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

    find_shedding_limits([tight_tuple, loose_tuple, dc_tuple, socwr_tuple, qcrm_tuple], filename, case_name; max_samples=max_samples, generator_order=ordered_gens, best_order=best_order)
    #group_plot(case_name, max_samples)
end

function group_plot(case_name, max_samples)
    saved_filename = "data/"*case_name
    plot_nflex_vs_y(saved_filename, "obj", "data/"*case_name*"_obj_samples"*string(max_samples); title=case_name*" obj vs flexible generators w/ max samples = "*string(max_samples))
    plot_nflex_vs_y(saved_filename, "obj", "data/"*case_name*"_obj_samples"*string(max_samples)*"_short"; title=case_name*" obj vs flexible generators w/ max samples = "*string(max_samples), tight_bounds=true)
    plot_nflex_vs_y(saved_filename, "total_pshed", "data/"*case_name*"_p_tot_samples"*string(max_samples); title=case_name*" total p shed vs flexible generators w/ max samples = "*string(max_samples))
    plot_nflex_vs_y(saved_filename, "max_pshed", "data/"*case_name*"_p_max_samples"*string(max_samples); title=case_name*" max p shed vs flexible generators w/ max samples = "*string(max_samples))
end

function plot_nflex_vs_y(
    csvfile::AbstractString,
    ycol::String,
    savepath::AbstractString;
    xlabel::AbstractString = "nflex",
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