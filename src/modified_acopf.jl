#Adds power generation fixing and load shedding




function build_augmented_opf(filename, flex, c_shed, pg_sol; backend=nothing)
    parsed_data = parse_matpower(filename; library=:pglib); 
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
        return ((p_shed=p_shed, q_shed=q_shed, p_shed_pos, p_shed_neg, q_shed_pos, q_shed_neg,), (c_pg_fix=c_pg_fix, c_p_shed_positive=c_p_shed_positive, c_p_shed_negative=c_p_shed_negative, c_q_shed_positive=c_q_shed_positive, c_q_shed_negative=c_q_shed_negative,))
    end
    return ac_opf_model(filename; user_callback = augment_opf, backend=backend)
end