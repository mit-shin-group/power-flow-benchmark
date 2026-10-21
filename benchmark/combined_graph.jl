using CSV
using DataFrames
using Plots
using Statistics
using OrderedCollections


function compute_ratios(file_list)
    ratios = OrderedDict{String, Vector{Float64}}()

    for file in file_list
        df = CSV.read(file, DataFrame)

        final_flex = df.flex_gen[end]
        groups = groupby(df, :initial_solution; sort=false)

        for g in groups
            idx = findfirst(x -> x < 1e-8, g.max_pshed)
            idx === nothing && continue

            init = String(g.initial_solution[1])
            push!(get!(ratios, init, Float64[]), g.flex_gen[idx] / final_flex)
        end
    end

    return ratios
end


function plot_cdfs(ratios, savepath; csvpath = savepath * ".csv", cost = false)
    plt = plot(legend = :bottomright)

    if cost
        xmax= 5

        x_label = "Percent Increase in Cost"
        y_label = "Portion of Cases at or Below Cost Increase"
    else
        xmax=1  
        x_label = "Portion of Generation Needed as Flexible"
        y_label = "Portion of Cases Feasible"
    end

    open(csvpath, "w") do io
        println(io, "Solution method,Avg")

        for (init, vals) in ratios
            sorted_vals = sort(vals)
            n = length(sorted_vals)

            # empirical CDF
            y = (1:n) ./ n

            # step plot
            plot!(plt, sorted_vals, y,
                  label = string(init),
                  seriestype = :steppost)

            # Average value for this solution
            println(io, "$(init),$(sum(vals) / length(vals))")
        end
    end

    xlabel!(x_label)
    ylabel!(y_label)
    ylims!(0, 1)
    xlims!(0, xmax)

    savefig(plt, savepath*".svg")
end


# ===== Usage =====

#ratios = compute_ratios(file_list)
#plt = plot_cdfs(ratios)
#display(plt)

##other options

#x axis - portion of slack generation available, y axis - average load shedding required


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
    grouped = groupby(all_data, :initial_solution; sort=false)

    plt = plot(yscale = :log10)

    for g in grouped
        # sort by flex_ratio for nicer curves
        g_sorted = sort(g, :flex_ratio)

        scatter!(
            plt,
            g_sorted.flex_ratio,
            g_sorted.pshed_ratio,
            label = first(g.initial_solution),
            markersize = 4
        )
    end

    xlabel!("flex_gen ratio")
    ylabel!("total_pshed ratio")
    title!("Flex Generation vs Load Shedding Tradeoff")
    ylims!(2e-6, 0.38)

    savefig(plt, savepath*".svg")
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

    savefig(plt, savepath * ".svg")
    return plt
end


#x axis - percent cost overage at feasible point, y axis - cumulative portion of cases at or below cost overage

function compute_ratios_cost_overage(file_list)
    ratios = OrderedDict{String, Vector{Float64}}()

    for file in file_list
        df = CSV.read(file, DataFrame)

        final_cost = df.obj[end]
        groups = groupby(df, :initial_solution; sort=false)  # order of appearance

        # Find the nflex value where each initial solution first becomes feasible
        feasibility_nflex = Float64[]

        for g in groups
            idx = findfirst(x -> x < 1.03e-8, g.max_pshed)

            if idx !== nothing
                push!(feasibility_nflex, g.nflex[idx])
            end
        end

        # Skip this file if no initial solution became feasible
        if isempty(feasibility_nflex)
            continue
        end

        # nflex at which the last initial solution becomes feasible
        worst_nflex = maximum(feasibility_nflex)

        # Evaluate every initial solution at worst_nflex
        for g in groups
            # Exact match first
            idx = findfirst(==(worst_nflex), g.nflex)

            # Otherwise use the first available nflex greater than worst_nflex
            if idx === nothing
                idx = findfirst(x -> x > worst_nflex, g.nflex)
            end

            # Skip this initial solution if it does not extend far enough
            idx === nothing && continue

            ratio = (g.obj[idx] - final_cost) / final_cost * 100

            

            init = String(g.initial_solution[1])

            if ratio > 10
                println(ratio, file, init)
            end
            push!(get!(ratios, init, Float64[]), ratio)
        end
    end

    return ratios
end


function compute_ratios_cost_overage_per_method(file_list)
    ratios = OrderedDict{String, Vector{Float64}}()

    for file in file_list
        df = CSV.read(file, DataFrame)

        final_cost = df.obj[end]
        groups = groupby(df, :initial_solution; sort=false)

        # For each solution method, find its own feasibility point and compute cost overage
        for g in groups
            # Find where this specific method first becomes feasible
            idx = findfirst(x -> x < 1.03e-8, g.max_pshed)

            # Skip if this method never becomes feasible
            if idx === nothing
                continue
            end

            # Compute cost overage at this method's own feasibility point
            ratio = (g.obj[idx] - final_cost) / final_cost * 100

            init = String(g.initial_solution[1])

            if ratio > 10
                println(ratio, " ", file, " ", init)
            end

            push!(get!(ratios, init, Float64[]), ratio)
        end
    end

    return ratios
end


function plot_flex_fraction_vs_generators(
    file_list::Vector{String},
    savepath
)
    # Store (total_generators, flex_fraction) for each initial solution
    data = Dict{Any, Tuple{Vector{Float64}, Vector{Float64}}}()

    for file in file_list
        df = CSV.read(file, DataFrame)

        # Total number of generators = final nflex entry
        total_generators = maximum(df.nflex)

        for g in groupby(df, :initial_solution)
            init = g.initial_solution[1]

            # First point with effectively zero load shedding
            idx = findfirst(x -> x < 1e-8, g.max_pshed)

            # Skip cases that never achieve feasibility
            if idx === nothing
                continue
            end

            nflex_needed = g.nflex[idx]
            flex_fraction = nflex_needed / total_generators

            if !haskey(data, init)
                data[init] = (Float64[], Float64[])
            end

            push!(data[init][1], total_generators)
            push!(data[init][2], flex_fraction)
        end
    end

    plt = plot()

    for (init, (generators, fractions)) in data
        # Sort by number of generators
        order = sortperm(generators)

        scatter!(
            plt,
            generators[order],
            fractions[order],
            marker = :circle,
            label = string(init)
        )
    end

    xlabel!("Total Number of Generators")
    ylabel!("Fraction of Flex Generation Needed")
    title!("Flex Generation Needed for Zero Load Shedding")
    ylims!(0, 1)

    savefig(plt, savepath * ".svg")

    return plt
end