using CSV
using DataFrames
using Plots
using ExaPowerIO

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

##other options

#x axis - portion of slack generation available, y axis - average load shedding required
#x axis - percent cost overage at feasible point, y axis - cumulative portion of cases at or below cost overage


using CSV, DataFrames, Plots

function plot_flex_vs_pshed(file_list::Vector{String}, loads::Vector{Float64}, savepath)
    @assert length(file_list) == length(loads) "file_list and loads must match"

    all_data = DataFrame()

    for (file, load) in zip(file_list, loads)
        df = CSV.read(file, DataFrame)

        # normalize flex_gen by max in that file
        max_flex = maximum(df.flex_gen)
        df.flex_ratio = df.flex_gen ./ max_flex

        # normalize total_pshed by provided load
        df.pshed_ratio = df.total_pshed ./ load

        append!(all_data, df)
    end

    # group by initial_solution
    grouped = groupby(all_data, :initial_solution)

    plt = plot()

    for g in grouped
        # sort by flex_ratio for nicer curves
        g_sorted = sort(g, :flex_ratio)

        plot!(
            plt,
            g_sorted.flex_ratio,
            g_sorted.pshed_ratio,
            label = first(g.initial_solution)
        )
    end

    xlabel!("flex_gen ratio")
    ylabel!("total_pshed ratio")
    title!("Flex Generation vs Load Shedding Tradeoff")

    savefig(plt, savepath*".png")
    return plt
end


function plot_flex_vs_pshed2(file_list::Vector{String}, loads::Vector{Float64}, savepath; interval::Float64=0.01)
    @assert length(file_list) == length(loads) "file_list and loads must match"

    all_data = DataFrame()

    for (file, load) in zip(file_list, loads)
        df = CSV.read(file, DataFrame)

        # normalize flex_gen by max in that file
        max_flex = maximum(df.flex_gen)
        df.flex_ratio = df.flex_gen ./ max_flex

        # normalize total_pshed by provided load
        df.pshed_ratio = df.total_pshed ./ load

        append!(all_data, df)
    end

    # create bins (0 to 1)
    all_data.flex_bin = floor.(all_data.flex_ratio ./ interval) .* interval

    # group by initial_solution and bin
    grouped = groupby(all_data, [:initial_solution, :flex_bin])

    # compute averages per bin
    agg = combine(grouped,
        :flex_ratio => mean => :flex_ratio_mean,
        :pshed_ratio => mean => :pshed_ratio_mean
    )

    # regroup only by initial_solution for plotting
    grouped_sol = groupby(agg, :initial_solution)

    plt = plot()

    for g in grouped_sol
        g_sorted = sort(g, :flex_ratio_mean)

        plot!(
            plt,
            g_sorted.flex_ratio_mean,
            g_sorted.pshed_ratio_mean,
            label = first(g.initial_solution)
        )
    end

    xlabel!("flex_gen ratio")
    ylabel!("total_pshed ratio")
    title!("Flex Generation vs Load Shedding Tradeoff")

    savefig(plt, savepath * ".png")
    return plt
end