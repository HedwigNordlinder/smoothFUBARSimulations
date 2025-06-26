using DataFrames, CSV, FASTX, Phylo, Plots, CodonMolecularEvolution, MolecularEvolution

function plot_tree(simulation_name)
    tree = read_newick_tree(simulation_name*".nwk")
    return plot(get_phylo_tree(tree))
end

function plot_loglog_rates(simulation_name)
    rate_frame = CSV.read(simulation_name*"_rates.csv", DataFrame)
    
    # Create scatter plot with log-transformed data
    p1 = scatter(log.(rate_frame.alphavec), log.(rate_frame.betavec),
                xlabel = "log(α)",
                ylabel = "log(β)", 
                title = "Log-Log Plot of Rates: "*simulation_name,
                label = "Data points",
                alpha = 0.7,
                markersize = 4)
    
    # Add y = x reference line
    x_range = xlims(p1)
    plot!(p1, [x_range[1], x_range[2]], [x_range[1], x_range[2]], 
          line = :dash, color = :red, linewidth = 2, label = "y = x")
    
    # Create histogram for alpha
    p2 = histogram(rate_frame.alphavec,
                   xlabel = "α",
                   ylabel = "Frequency",
                   title = "Distribution of α: "*simulation_name,
                   alpha = 0.7,
                   color = :blue,
                   bins = 30)
    
    # Create histogram for beta
    p3 = histogram(rate_frame.betavec,
                   xlabel = "β",
                   ylabel = "Frequency", 
                   title = "Distribution of β: "*simulation_name,
                   alpha = 0.7,
                   color = :green,
                   bins = 30)
    
    return p1, p2, p3
end

function plot_scenario_information(simulation_name)
    scenario = load_scenario(simulation_name * "_scenario.json")
    Ne = t -> effective_population_size(scenario, t)
    s = t -> sampling_rate(scenario, t)  # renamed to avoid conflict
    
    # Create time vector from 0 to 1
    t_vec = range(0, 100, length=1000)
    
    # Evaluate functions over time
    Ne_values = [Ne(t) for t in t_vec]
    sampling_values = [s(t) for t in t_vec]
    
    # Create first plot - Effective Population Size
    p1 = plot(t_vec, Ne_values,
              xlabel="Time",
              ylabel="Effective Population Size (Ne)",
              title="Effective Population Size Over Time",
              linewidth=2,
              color=:blue,
              grid=true)
    
    # Create second plot - Sampling Rate
    p2 = plot(t_vec, sampling_values,
              xlabel="Time", 
              ylabel="Sampling Rate",
              title="Sampling Rate Over Time",
              linewidth=2,
              color=:red,
              grid=true)
    
    return p1, p2
end

function save_tree_report(simulation_name, output_filename=nothing)
    # Set default filename if not provided
    if output_filename === nothing
        output_filename = simulation_name * "_tree_report.pdf"
    end
    
    # Generate all plots
    tree_plot = plot_tree(simulation_name)
    scatter_plot, alpha_hist, beta_hist = plot_loglog_rates(simulation_name)
    scenario_p1, scenario_p2 = plot_scenario_information(simulation_name)
    
    # Create a combined layout with all plots
    # Tree on top row
    # Three rate plots in middle row 
    # Two scenario plots in bottom row
    combined_plot = plot(tree_plot, 
                        scatter_plot, alpha_hist, beta_hist,
                        scenario_p1, scenario_p2,
                        layout = @layout([a; [b c d]; [e f]]),
                        size = (1200, 1200),
                        plot_title = "Tree Report: " * simulation_name,
                        plot_titlefontsize = 16)
    
    # Save to PDF
    savefig(combined_plot, output_filename)
    
    println("Tree report saved to: ", output_filename)
    return output_filename
end