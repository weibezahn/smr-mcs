##### plots #####
# requires results from the main file

using CairoMakie
include("functions_plots.jl")

##### comparison plot for Roulstone vs. Rothwell scaling #####

# α range
x = 0:0.01:1;
# β Roulstone
y = x;
# β Rothwell
z = 1 .+ log.(x)/log(2);

# define plot
fig_theory = Figure();

ax_theory = Axis(fig_theory[1,1], xlabel = "α", ylabel = "β(α)");
xlims!(0, 1)
ylims!(-2, 1)

roulstone = lines!(x, y, label = "Roulstone", linewidth = 3, color = :darkblue);
rothwell = lines!(x, z, lable = "Rothwell", linewidth = 3, color = :green);
hlines!(ax_theory, [0], color = :gray);

Legend(fig_theory[1, 1],
    [roulstone, rothwell],
    [L"β^\text{Roulstone}", L"β^\text{Rothwell}"],
    tellheight = false,
    tellwidth = false,
    halign = :right, valign = :bottom,
    framevisible = false, orientation = :vertical);

fig_theory
save("$outputpath/fig-theory.pdf", fig_theory);

##### comparison plot for investment cost from manufacturers vs. estimation #####

# choose scaling parameters for the plot
scaling_plot = [0.20, 0.75];

fig_invest_comparison = investment_plot(pjs, scaling_plot)
save("$outputpath/fig-investment_comparison.pdf", fig_invest_comparison);

##### histogram plots for comparison of estimation approaches #####

hist_invest = Figure();

for i in 1:3, j in 1:5
    hist_invest_plot(n, wacc, electricity_price, pjs[j+5*(i-1)], i, j, hist_invest)
end

Legend(hist_invest[4,1:5],
    [roulstone, rothwell],
    ["Roulstone", "Rothwell"],
    framevisible = false, orientation = :horizontal)

hist_invest
save("$outputpath/fig-histogram_investment.pdf", hist_invest);

##### probability density plots for comparison of estimation approaches #####

density_invest = Figure();

for i in 1:3, j in 1:5
    density_invest_plot(n, wacc, electricity_price, pjs[j+5*(i-1)], i, j, density_invest)
end

Legend(density_invest[4,1:5],
    [roulstone, rothwell],
    ["Roulstone", "Rothwell"],
    framevisible = false, orientation = :horizontal)

density_invest
save("$outputpath/fig-density_investment.pdf", density_invest);

##### boxplots Monte Carlo simulation results #####
# requires results for all 15 reactor concepts

fig_mcs_npv = mcs_plot(npv_results, "NPV", "[USD/MW]")
fig_mcs_lcoe = mcs_plot(lcoe_results, "LCOE", "[USD/MWh]")

save("$outputpath/fig-mcs_npv-$opt_scaling.pdf", fig_mcs_npv);
save("$outputpath/fig-mcs_lcoe-$opt_scaling.pdf", fig_mcs_lcoe);

##### heatmaps sensitivity indices #####
# requires sensitvity results

fig_si_npv = si_plot(si_npv_results, "NPV")
fig_si_lcoe = si_plot(si_lcoe_results, "LCOE")

save("$outputpath/fig-si_npv-$opt_scaling.pdf", fig_si_npv);
save("$outputpath/fig-si_lcoe-$opt_scaling.pdf", fig_si_lcoe);

##### lcoe comparison plot #####
# requires results for all 15 reactor concepts

using CSV, DataFrames

# read LCOE data from CSV into a dataframe
lcoe_dat = CSV.File("$inputpath/lcoe_data.csv") |> DataFrame;

# read LCOE data from project simulations
lcoe_bounds = select(lcoe_summary, [:q25, :q75]);

# collect plot data
lcoe_plot_data = vcat(
    select(lcoe_dat, [:technology, :lower_bound, :upper_bound]),
    DataFrame(technology = ["BWR & PWR SMRs", "HTR SMRs", "SFR SMRs"],
    lower_bound = [minimum(lcoe_bounds[1:9,:q25]), minimum(lcoe_bounds[10:12,:q25]), minimum(lcoe_bounds[13:15,:q25])],
    upper_bound = [maximum(lcoe_bounds[1:9,:q75]), maximum(lcoe_bounds[10:12,:q75]), maximum(lcoe_bounds[13:15,:q75])])
);

# define LCOE plot
if opt_scaling == "manufacturer"
    plot_scaling = "Manufacturer";
elseif opt_scaling == "roulstone"
    plot_scaling = "Roulstone";
elseif opt_scaling == "rothwell"
    plot_scaling = "Rothwell";
elseif opt_scaling == "uniform"
    plot_scaling = "uniform";
else
    @error("scaling not defined")
end
xlabel = "[USD/MWh]";
yticks = lcoe_plot_data[!,:technology];

col = vcat(fill(1,8),fill(2,4),3,4,5);
colormap = [:darkgreen, :darkblue];

fig_lcoe_comparison = Figure();
ax_lcoe = Axis(fig_lcoe_comparison[1,1], yticks = (1:length(yticks), yticks), xscale = log10, xlabel = xlabel);

xlims!(10, 25000)

rangebars!(ax_lcoe, 1:length(yticks), lcoe_plot_data[!,2], lcoe_plot_data[!,3], linewidth = 6, whiskerwidth = 8, direction = :x, color = col);
hlines!(ax_lcoe, [8.5, 12.5], linestyle = :dash, color = :red);
text!([15000,15000,16], [4, 10.5, 14]; text = ["Renewables\n(LAZARD)", "Conventionals\n(LAZARD)", "SMR Tech.\n($plot_scaling)"], align = (:center, :center), justification = :center, rotation = π/2);

text!(lcoe_plot_data[!,2], 1:length(yticks), text = string.(round.(Int,lcoe_plot_data[!,2])), align = (:right, :center), offset = (-10,0));
text!(lcoe_plot_data[!,3], 1:length(yticks), text = string.(round.(Int,lcoe_plot_data[!,3])), align = (:left, :center), offset = (10,0));

Label(fig_lcoe_comparison[1, 1, Top()], "LCOE Comparison", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

fig_lcoe_comparison
save("$outputpath/fig-lcoe_comparison-$opt_scaling.pdf", fig_lcoe_comparison);