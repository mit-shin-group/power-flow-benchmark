using CSV
using DataFrames
using Plots

function compute_ratios(file_list)
    # Store ratios grouped by initial_solution
    ratios = Dict{Any, Vector{Float64}}()

    for file in file_list
        df = CSV.read(file, DataFrame)

        # final flex_gen (last row of CSV)
        final_flex = df.flex_gen[end]

        # group by initial_solution
        for g in groupby(df, "initial_solution")
            init = g.initial_solution[1]

            # find first index where total_pshed < 1e-8
            idx = findfirst(x -> x < 1e-8, g.max_pshed)

            if idx !== nothing
                ratio = g.flex_gen[idx] / final_flex

                if !haskey(ratios, init)
                    ratios[init] = Float64[]
                end
                push!(ratios[init], ratio)
            end
        end
    end

    return ratios
end


function plot_cdfs(ratios, savepath)
    plt = plot()

    for (init, vals) in ratios
        sorted_vals = sort(vals)
        n = length(sorted_vals)

        # empirical CDF
        y = (1:n) ./ n

        # step plot (piecewise constant)
        plot!(plt, sorted_vals, y,
              label = string(init),
              seriestype = :steppost)
    end

    xlabel!("Portion of Generation Needed as Slack")
    ylabel!("Portion of Cases Feasible")
    xlims!(0, 1)
    ylims!(0, 1)

    savefig(plt, savepath*".png")
end


# ===== Usage =====

#ratios = compute_ratios(file_list)
#plt = plot_cdfs(ratios)
#display(plt)