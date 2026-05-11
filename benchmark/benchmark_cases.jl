using ExaPowerIO

short_cases = [
"pglib_opf_case3_lmbd", 
"pglib_opf_case5_pjm",
"pglib_opf_case14_ieee",
"pglib_opf_case24_ieee_rts",
"pglib_opf_case30_as",]

full_cases = [
"pglib_opf_case3_lmbd", 
"pglib_opf_case5_pjm",
"pglib_opf_case14_ieee",
"pglib_opf_case24_ieee_rts",
"pglib_opf_case30_as",
"pglib_opf_case30_ieee",
"pglib_opf_case39_epri",
"pglib_opf_case57_ieee",
"pglib_opf_case60_c",
"pglib_opf_case73_ieee_rts",
"pglib_opf_case89_pegase",
"pglib_opf_case118_ieee",
"pglib_opf_case162_ieee_dtc",
"pglib_opf_case179_goc",
"pglib_opf_case200_activ",
"pglib_opf_case240_pserc",
"pglib_opf_case300_ieee",
"pglib_opf_case500_goc",
"pglib_opf_case588_sdet",
"pglib_opf_case793_goc",
"pglib_opf_case1354_pegase",
"pglib_opf_case1803_snem",
"pglib_opf_case1951_rte",
"pglib_opf_case2000_goc",
"pglib_opf_case2312_goc",
"pglib_opf_case2383wp_k",
"pglib_opf_case2736sp_k",
"pglib_opf_case2737sop_k",
"pglib_opf_case2742_goc",
"pglib_opf_case2746wop_k",
"pglib_opf_case2746wp_k",
"pglib_opf_case2848_rte",
"pglib_opf_case2853_sdet",
"pglib_opf_case2868_rte",
"pglib_opf_case2869_pegase",
"pglib_opf_case3012wp_k",
"pglib_opf_case3022_goc",
"pglib_opf_case3120sp_k",
"pglib_opf_case3375wp_k",
"pglib_opf_case3970_goc",
"pglib_opf_case4020_goc",
"pglib_opf_case4601_goc",
"pglib_opf_case4619_goc",
"pglib_opf_case4661_sdet",
"pglib_opf_case4837_goc",
"pglib_opf_case4917_goc",
"pglib_opf_case5658_epigrids",
"pglib_opf_case6468_rte",
"pglib_opf_case6470_rte",
"pglib_opf_case6495_rte",
"pglib_opf_case6515_rte",
"pglib_opf_case7336_epigrids",
"pglib_opf_case8387_pegase",
"pglib_opf_case9241_pegase",
"pglib_opf_case9591_goc",
"pglib_opf_case10000_goc",
"pglib_opf_case10192_epigrids",
"pglib_opf_case10480_goc",
"pglib_opf_case13659_pegase",
"pglib_opf_case19402_goc",
"pglib_opf_case20758_epigrids",
"pglib_opf_case24464_goc",
"pglib_opf_case30000_goc",
"pglib_opf_case78484_epigrids",]

function benchmark_ac_vs_dc(cases; max_samples = 100, sort_order = :best, form = :standard)
    for case in cases
        if form == :api
            case = case*"__api"
            save_folder = "saved_data/api/"
        else
            save_folder = "saved_data/"
        end
        case_label = split(case, "pglib_opf_")[2]
        if sort_order == :best
            case_label_saved = case_label*"_best_ordered"
        end
        case_label_saved = case_label_saved*"_samples"*string(max_samples)

        println(case_label_saved)
        println(case_label)
        println(case)
        if !isfile(save_folder*case_label_saved)
            test_ac_dc_combo("data/"*case*".m", case_label; max_samples = max_samples, sort_order = sort_order, save_folder=save_folder)
        end
    end
end

function make_plots(savefile, savefile2, savefile3; cases=full_cases, max_samples=100, sort_order=:best, interval = 0.01)
    case_list = String[]
    load_list = Float64[]
    for case in cases
        total_load = sum(b.pd for b in ExaPowerIO.parse_matpower("data/"*case*".m";).bus)
        case_label = split(case, "pglib_opf_")[2]
        if sort_order == :best
            case_label = case_label* "_best_ordered"
        end
        case_label = case_label*"_samples"*string(max_samples)
        if isfile("saved_data/"*case_label)
            push!(case_list, "saved_data/"*case_label)
            push!(load_list, total_load)
        end
    end
    ratios = compute_ratios(case_list)
    #plt = plot_cdfs(ratios, "saved_data/formatted_data/"*savefile)
    cost_ratios = compute_ratios_cost_overage(case_list)
    #plt = plot_cdfs(cost_ratios, "saved_data/formatted_data/"*savefile3; cost = true)
    #plot_flex_vs_pshed2(case_list, load_list, "saved_data/formatted_data/"*savefile2; interval=interval)
    plot_flex_vs_pshed(case_list, load_list, "saved_data/formatted_data/"*savefile2;)
end