using DataFrames, CSV, Plots, ExaModels, ExaModelsPower, MadNLP, MadNLPGPU, ExaPowerIO, CUDA

function find_shedding_limits(pg_solutions, filename, save_filename; max_samples=30, gen_range= [])
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

    generator_order = reverse(sort(1:length(gens), by = i -> gens[i].pmax))

    flex_generators = Int[]


    if isempty(gen_range)
        if L_gen ≤ max_samples
            gen_range = collect(0:L_gen+1)
        else
            gen_range = unique(round.(Int, range(0, L_gen+1, length=max_samples)))
        end
    end
    
    for sol_tuple in pg_solutions
        flex_generators = Int[]
        for gen in gen_range
            j=0
            while j < 3
            
                println(length(flex_generators))
                println(sol_tuple.initial_method)
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
            if isempty(flex_generators)
                flex_generators = [1]
            else
                append!(flex_generators, (maximum(flex_generators) + 1):gen)
            end
        end
    end

    CSV.write("data/"*save_filename, df)
end

function test_ac_dc_combo(filename, case_name; max_samples=30, ac_order = false)
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


    saved_filename = "data/"*case_name

    if ac_order
        ac_ordered_gens = find_ac_gen_order(tight_tuple.pg_sol)
        case_name = case_name*"_ac_ordered"
    else
        ac_ordered_gens = []
    end

    find_shedding_limits([tight_tuple, loose_tuple, dc_tuple, socwr_tuple, qcrm_tuple], filename, case_name; max_samples=max_samples, gen_range=ac_ordered_gens)
    plot_nflex_vs_y(saved_filename, "obj", "data/"*case_name*"_obj_samples"*string(max_samples); title=case_name*" obj vs flexible generators w/ max samples = "*string(max_samples))
    plot_nflex_vs_y(saved_filename, "total_pshed", "data/"*case_name*"_p_tot_samples"*string(max_samples); title=case_name*" total p shed vs flexible generators w/ max samples = "*string(max_samples))
    plot_nflex_vs_y(saved_filename, "max_pshed", "data/"*case_name*"_p_max_samples"*string(max_samples); title=case_name*" max p shed vs flexible generators w/ max samples = "*string(max_samples))
    println(ac_ordered_gens)
    println(tight_tuple.pg_sol)
end

function plot_nflex_vs_y(
    csvfile::AbstractString,
    ycol::String,
    savepath::AbstractString;
    xlabel::AbstractString = "nflex",
    ylabel::AbstractString = ycol,
    title::AbstractString = ""
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

    for (i, sol) in enumerate(unique(df.initial_solution))
        sub = df[df[!, "initial_solution"] .== sol, :]
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