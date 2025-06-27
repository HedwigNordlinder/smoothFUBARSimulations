using Pkg
Pkg.activate(".")
Pkg.add(url="https://github.com/HedwigNordlinder/CodonMolecularEvolution.jl")
using MolecularEvolution, CodonMolecularEvolution, Distributions, Plots, Phylo, DataFrames, CSV
include("plot_simulation_data.jl")

α_distribution = Gamma(10,0.1)
β_distribution = Exponential(1)
function create_simulation_parameter_csv(output_file::String = "simulation_parameters.csv")
    # Parameter combinations to explore
    ntaxa_values = [10, 20, 50, 100, 200]
    nsites_values = [100, 500, 1000]
    diversifying_sites_values = [1, 10, 50]
    scenarios = ["LogisticScenario()", "SeasonalScenario(;sin_divisor=1.0)"]
    
    # Static parameters (adjust these based on your uploaded CSV)
    static_params = (
        alpha_distribution = "Gamma(10,0.1)",
        beta_distribution = "Exponential(1)", 
        nucleotide_model = "default",
        f3x4_model = "default",
        target_normalisation = 10.0
    )
    
    # Create all combinations
    combinations = collect(Iterators.product(ntaxa_values, nsites_values, diversifying_sites_values, scenarios))
    
    # Create DataFrame
    df = DataFrame(
        scenario_name = String[],
        ntaxa = Int[],
        nsites = Int[],
        diversifying_sites = Int[],
        coalescence_scenario = String[],
        alpha_distribution = String[],
        beta_distribution = String[],
        nucleotide_model = String[],
        f3x4_model = String[],
        target_normalisation = Float64[]
    )
    
    # Generate rows
    for (ntaxa, nsites, div_sites, scenario) in combinations
        scenario_name = "$(lowercase(replace(scenario, "()" => "")))_t$(ntaxa)_s$(nsites)_d$(div_sites)"
        
        push!(df, (
            scenario_name = scenario_name,
            ntaxa = ntaxa,
            nsites = nsites,
            diversifying_sites = div_sites,
            coalescence_scenario = scenario,
            alpha_distribution = static_params.alpha_distribution,
            beta_distribution = static_params.beta_distribution,
            nucleotide_model = static_params.nucleotide_model,
            f3x4_model = static_params.f3x4_model,
            target_normalisation = static_params.target_normalisation
        ))
    end
    
    # Write to CSV
    CSV.write(output_file, df)
    
    println("Created $(nrow(df)) simulation scenarios in $output_file")
    println("Combinations: $(length(ntaxa_values)) ntaxa × $(length(nsites_values)) nsites × $(length(diversifying_sites_values)) diversifying_sites × $(length(scenarios)) scenarios")
    
    return df
end

# Run it
df = create_simulation_parameter_csv("simulation_parameters.csv")
# This is the main simulation function

#result_star = simulate_k_diversifying_sites(ntaxa,SeasonalScenario(;sampling_divisor=25), α_distribution, β_distribution, 
#                                        nsites, diversifying_sites, CodonMolecularEvolution.demo_nucmat, 
#                                        CodonMolecularEvolution.demo_f3x4)

#save_simulation_data(result_star, name="seasonal")     
#save_tree_report("seasonal")                               
#result_ladder = simulate_k_diversifying_sites(ntaxa,LogisticScenario(), α_distribution, β_distribution, 
#                                        nsites, diversifying_sites, CodonMolecularEvolution.demo_nucmat, 
#                                        CodonMolecularEvolution.demo_f3x4)

#save_simulation_data(result_ladder, name="logistic")   

 
#save_tree_report("logistic")
run_simulation_batch("simulation_parameters.csv","simulations")