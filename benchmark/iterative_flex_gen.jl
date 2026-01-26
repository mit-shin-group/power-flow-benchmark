using DataFrames, CSV, Plots, ExaModels, ExaModelsPower, MadNLP, MadNLPGPU, ExaPowerIO, CUDA

function find_shedding_limits(pg_solutions, filename, save_filename; max_samples=30)
    parsed_data = parse_matpower(filename; library=:pglib)
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


    if L_gen ≤ max_samples
        gen_range = collect(0:L_gen)
    else
        gen_range = unique(round.(Int, range(0, L_gen, length=max_samples)))
        
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
                    push!(df, [sol_tuple.initial_method, length(flex_generators), result.objective, sum(p_sheds), maximum(p_sheds), sum(q_sheds), maximum(q_sheds)])
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

function test_ac_dc_combo(filename, case_name; max_samples=30)
    m_ac, v_ac, c_ac = ac_opf_model(filename; backend = CUDABackend())
    result_loose = madnlp(m_ac, tol=1e-4)
    result_tight = madnlp(m_ac, tol=1e-8)
    m_dc, v_dc, c_dc = dcopf_model(filename; backend=CUDABackend())
    result_dc = madnlp(m_dc, tol=1e-8)
    loose_tuple = (pg_sol=Array(solution(result_loose, v_ac.pg)), initial_method = "ACOPF tol 1e-4")
    tight_tuple = (pg_sol=Array(solution(result_tight, v_ac.pg)), initial_method = "ACOPF tol 1e-8")
    dc_tuple = (pg_sol=Array(solution(result_dc, v_dc.pg)), initial_method = "DCOPF tol 1e-8")
    saved_filename = "data/"*case_name
    find_shedding_limits([tight_tuple, loose_tuple, dc_tuple], filename, case_name; max_samples=max_samples)
    plot_nflex_vs_y(saved_filename, "obj", "data/"*case_name*"_obj_samples"*string(max_samples); title=case_name*" obj vs flexible generators w/ max samples = "*string(max_samples))
    plot_nflex_vs_y(saved_filename, "total_pshed", "data/"*case_name*"_p_tot_samples"*string(max_samples); title=case_name*" total p shed vs flexible generators w/ max samples = "*string(max_samples))
    plot_nflex_vs_y(saved_filename, "max_pshed", "data/"*case_name*"_p_max_samples"*string(max_samples); title=case_name*" max p shed vs flexible generators w/ max samples = "*string(max_samples))
end

function plot_nflex_vs_y(
    csvfile::AbstractString,
    ycol::String,
    savepath::AbstractString;
    xlabel::AbstractString = "nflex",
    ylabel::AbstractString = ycol,
    title::AbstractString = ""
)

    df = CSV.read(csvfile, DataFrame)


    @assert "initial_solution" in names(df)
    @assert "nflex" in names(df)
    @assert ycol in names(df)

    xmin = minimum(df.nflex)
    xmax = maximum(df.nflex)

    # ensure strictly positive width
    if xmin == xmax
        xmax = xmin + 1
    end


    p = plot(
        xlabel = xlabel,
        ylabel = ylabel,
        yscale = :log10,
        title = title,
        legend = :best,
        #xlims = (xmin, xmax),
    )

    for sol in unique(df.initial_solution)
        sub = df[df[!,"initial_solution"] .== sol, :]
        sort!(sub, "nflex")

    
        plot!(
            p,
            sub[!, "nflex"],
            sub[!, ycol],
            label = sol
        )
    end

    savefig(p, savepath)
    return p
end