cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using PESGM_LineDynamics

const ETL = PESGM_LineDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("9bus_analysis_functions.jl")

function retrieve_samples(filename)
    samples_cut = DataFrame(CSV.File(filename))
    return samples_cut
end

function generate_new_samples(N, seed, exp_folder)
    samples = generate_single_case_gain_params() # generate all combinations 
    samples_cut = get_random_sample(samples, N, seed) # randomly sample N 
    mkdir(exp_folder)
    CSV.write(exp_folder*"/samples.csv", samples_cut)
    return samples_cut
end

function generate_data(exp_df, samples_cut, exp_folder)
    for i = 1:size(exp_df)[1];
        line_scale = exp_df[i,:line_scale]
        inv_share = exp_df[i,:inv_share]
        load_scale = exp_df[i,:load_scale]
        gfl_share = exp_df[i,:gfl_share]
        case = exp_df[i,:case_name]

        generate_single_case_3(samples_cut, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, exp_folder, case)

    end
end

function generate_figures(exp_df, gain_case, exp_folder, save_folder)
    # Define yaxis scales for each case 
    case1_scale = 0.1;
    case2_scale = 0.02;
    case3_scale = 0.4;
    case4_scale = 0.2;

    yscales = [case1_scale, case1_scale, case2_scale, case2_scale*10, case3_scale, case3_scale, case4_scale, case4_scale]
    legend_key = [:bottomright, false, false, false, false, false, false, false]
    legend_key = [false, false, :bottomright, false, false, false, false, false]

    ylabel_key = [L"\Re \ (\lambda)", L"\Re \ (\lambda)", "", "", "", "", "", ""]

    # Make boxplots for each case 
    for i = 1:size(exp_df)[1];
        load_scale = exp_df[i,:load_scale]
        case = exp_df[i,:case_name]
        case_name = case*"_"*string(load_scale)
        suffix = case*"_"*string(load_scale)*".csv"
        alg = DataFrame(CSV.File(exp_folder*"/alg_"*suffix))
        dyn = DataFrame(CSV.File(exp_folder*"/dyn_"*suffix))
        mssb = DataFrame(CSV.File(exp_folder*"/mssb_"*suffix))

        make_case_boxplot_3(alg, dyn, mssb, legend_key[i], ylabel_key[i], case_name, save_folder, yscales[i])
    end

    # choose case to make gain plots for 
    i = gain_case 
    load_scale = exp_df[i,:load_scale]
    case = exp_df[i,:case_name]

    case_name = case*"_"*string(load_scale)
    suffix = case*"_"*string(load_scale)*".csv"
    alg = DataFrame(CSV.File(exp_folder*"/alg_"*suffix))
    dyn = DataFrame(CSV.File(exp_folder*"/dyn_"*suffix))
    mssb = DataFrame(CSV.File(exp_folder*"/mssb_"*suffix))

    make_gain_boxplots([alg, dyn, mssb], [L"\mathrm{statpi}", L"\mathrm{dynpi}", L"\mathrm{MSSB}"], :inside, case_name, save_folder)
end


file_name = "../data/json_data/9bus_VSM3_SM1_GFL2.json"
#N = 1000; # Number of samples of parameters 
#R = 1234; # Random seed for doing selection of parameters 
exp_folder = "Final" # name of folder for this experiment - data is saved here
save_folder = "Revised_figs" # figures are saved here 

# Load df of cases 
exp_df = DataFrame(CSV.File("pesgm_exp_df.csv"))
samples_cut = retrieve_samples("Final/samples.csv")

#generate_data(exp_df, samples_cut, exp_folder)

generate_figures(exp_df, 1, exp_folder, save_folder)
