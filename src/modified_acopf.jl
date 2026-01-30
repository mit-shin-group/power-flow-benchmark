#Adds power generation fixing and load shedding

using PowerModels, Ipopt


function augmented_opf_L1(filename, c_shed, pg_sol, lambda; backend=nothing)
    parsed_data = ExaPowerIO.parse_matpower(filename;)# library=:pglib); 
    L_pg = length(parsed_data.gen)
    L_bus = length(parsed_data.bus)

    pg_sol = Array(pg_sol)

    bus_list = [(i=i,) for i in 1:L_bus]

    function augment_opf(core, vars, cons)
        p_shed = variable(core, L_bus;)
        p_shed_pos = variable(core, L_bus; lvar = 0)
        p_shed_neg = variable(core, L_bus; lvar = 0)
        q_shed = variable(core, L_bus;)
        q_shed_pos = variable(core, L_bus; lvar = 0)
        q_shed_neg = variable(core, L_bus; lvar = 0)
        o_shed = objective(core, c_shed*(p_shed_pos[i]+q_shed_pos[i]+p_shed_neg[i]+q_shed_neg[i]) for i in 1:L_bus)
        c_active_power_balance_shed = constraint!(core, cons.c_active_power_balance, b.i => -p_shed[b.i] for b in bus_list)
        c_reactive_power_balance_shed = constraint!(core, cons.c_reactive_power_balance, i => -q_shed[i] for i in 1:L_bus)

        c_p_shed_positive = constraint(core, p_shed_pos[i]-p_shed[i] for i in 1:L_bus; ucon=fill(Inf, L_bus))
        c_p_shed_negative = constraint(core, -p_shed_neg[i]-p_shed[i] for i in 1:L_bus; lcon=fill(-Inf, L_bus))
        c_q_shed_positive = constraint(core, q_shed_pos[i]-q_shed[i] for i in 1:L_bus; ucon=fill(Inf, L_bus))
        c_q_shed_negative = constraint(core, -q_shed_neg[i]-q_shed[i] for i in 1:L_bus; lcon=fill(-Inf, L_bus))

        pg = vars.pg

        #Add L1 regression
        p_diff = variable(core, L_pg)
        p_diff_pos = variable(core, L_pg; lvar=0)
        p_diff_neg = variable(core, L_pg; lvar=0)

        println(pg_sol)

        fixed_generators = [
                            (i = i,
                            g = pg_sol[i])
                            for i in 1:L_pg
                        ]
        c_p_diff = constraint(core, pg[fg.i] - fg.g - p_diff[fg.i] for fg in fixed_generators)
        c_p_diff_positive = constraint(core, p_diff_pos[i]-p_diff[i] for i in 1:L_pg; ucon=fill(Inf, L_pg))
        c_p_diff_negative = constraint(core, -p_diff_neg[i]-p_diff[i] for i in 1:L_pg; lcon=fill(-Inf, L_pg))

        o_p_diff = objective(core, lambda*(p_diff_neg[i]+p_diff_pos[i]) for i in 1:L_pg)

        return ((p_shed=p_shed, q_shed=q_shed, p_shed_pos=p_shed_pos, p_shed_neg=p_shed_neg, q_shed_pos=q_shed_pos, q_shed_neg=q_shed_neg, p_diff=p_diff, p_diff_pos=p_diff_pos, p_diff_neg=p_diff_neg,), 
        (c_p_shed_positive=c_p_shed_positive, c_p_shed_negative=c_p_shed_negative, c_q_shed_positive=c_q_shed_positive, c_q_shed_negative=c_q_shed_negative, c_p_diff_positive=c_p_diff_positive, c_p_diff_negative=c_p_diff_negative,)) #c_p_diff=c_p_diff, 
    end
    return ac_opf_model(filename; user_callback = augment_opf, backend=backend)
end


function build_augmented_opf(filename, flex, c_shed, pg_sol; backend=nothing)
    parsed_data = ExaPowerIO.parse_matpower(filename;)# library=:pglib); 
    L_pg = length(parsed_data.gen)
    L_bus = length(parsed_data.bus)

    pg_sol = Array(pg_sol)

    fixed_generators = [
                            (i = i,
                            g = clamp(pg_sol[i],
                                    parsed_data.gen[i].pmin,
                                    parsed_data.gen[i].pmax))
                            for i in setdiff(1:L_pg, flex)
                        ]


    bus_list = [(i=i,) for i in 1:L_bus]

    function augment_opf(core, vars, cons)
        p_shed = variable(core, L_bus;)
        p_shed_pos = variable(core, L_bus; lvar = 0)
        p_shed_neg = variable(core, L_bus; lvar = 0)
        q_shed = variable(core, L_bus;)
        q_shed_pos = variable(core, L_bus; lvar = 0)
        q_shed_neg = variable(core, L_bus; lvar = 0)
        o_shed = objective(core, c_shed*(p_shed_pos[i]+q_shed_pos[i]+p_shed_neg[i]+q_shed_neg[i]) for i in 1:L_bus)
        c_active_power_balance_shed = constraint!(core, cons.c_active_power_balance, b.i => -p_shed[b.i] for b in bus_list)
        c_reactive_power_balance_shed = constraint!(core, cons.c_reactive_power_balance, i => -q_shed[i] for i in 1:L_bus)

        c_p_shed_positive = constraint(core, p_shed_pos[i]-p_shed[i] for i in 1:L_bus; ucon=fill(Inf, L_bus))
        c_p_shed_negative = constraint(core, -p_shed_neg[i]-p_shed[i] for i in 1:L_bus; lcon=fill(-Inf, L_bus))
        c_q_shed_positive = constraint(core, q_shed_pos[i]-q_shed[i] for i in 1:L_bus; ucon=fill(Inf, L_bus))
        c_q_shed_negative = constraint(core, -q_shed_neg[i]-q_shed[i] for i in 1:L_bus; lcon=fill(-Inf, L_bus))

        pg = vars.pg
        c_pg_fix = constraint(core, pg[fg.i] - fg.g for fg in fixed_generators)
        return ((p_shed=p_shed, q_shed=q_shed, p_shed_pos=p_shed_pos, p_shed_neg=p_shed_neg, q_shed_pos=q_shed_pos, q_shed_neg=q_shed_neg,), (c_pg_fix=c_pg_fix, c_p_shed_positive=c_p_shed_positive, c_p_shed_negative=c_p_shed_negative, c_q_shed_positive=c_q_shed_positive, c_q_shed_negative=c_q_shed_negative,))
    end
    return ac_opf_model(filename; user_callback = augment_opf, backend=backend)
end

function powermodels_pg_result(filename, model_type)
    result = solve_opf(filename, model_type, Ipopt.Optimizer)
    gen_results = Float64[]
    for i in sort(parse.(Int, keys(result["solution"]["gen"])))
        append!(gen_results, result["solution"]["gen"][string(i)]["pg"])
    end
    return gen_results
end