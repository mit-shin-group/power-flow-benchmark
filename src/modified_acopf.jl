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

    function augment_opf(core, vars, cons)
        p_shed = variable(core, L_bus;)
        q_shed = variable(core, L_bus;)
        o_shed = objective(core, c_shed*(p_shed[i]^2+q_shed[i]^2) for i in 1:L_bus)
        c_active_power_balance_shed = constraint!(core, cons.c_active_power_balance, i => -p_shed[i] for i in 1:L_bus)
        c_reactive_power_balance_shed = constraint!(core, cons.c_reactive_power_balance, i => -q_shed[i] for i in 1:L_bus)

        pg = vars.pg
        c_pg_fix = constraint(core, pg[fg.i] - fg.g for fg in fixed_generators)
        return ((p_shed=p_shed, q_shed=q_shed,), (c_pg_fix=c_pg_fix,))
    end
    return ac_opf_model(filename; user_callback = augment_opf, backend=backend)
end