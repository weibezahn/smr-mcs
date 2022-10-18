##### plots #####

using StatsPlots

##### NPV plot #####
boxplot(Matrix(npv_results), legend=false,
    xticks=(1:15, names(npv_results)), xrotation = 45,
    ylabel="NPV [USD/MW]",
    colour=:white)
vline!([1.5,9.5,12.5,15.5], colour=:red)
annotate!(1.25, -7e8, text("BWR type SMR concepts",6, rotation=90))
annotate!(9.25, -7e8, text("PWR type SMR concepts",6, rotation=90))
annotate!(12.25, -7e8, text("HTR type SMR concepts",6, rotation=90))
annotate!(15.25, -7e8, text("SFR type SMR concepts",6, rotation=90))
savefig("fig-npv_mcs.pdf")

##### LCOE plot #####
boxplot(Matrix(lcoe_results), legend=false,
    xticks=(1:15, names(lcoe_results)), xrotation = 45,
    ylabel="LCOE [USD/MWh]",
    colour=:white)
vline!([1.5,9.5,12.5,15.5], colour=:red)
annotate!(1.25, 1.25e4, text("BWR type SMR concepts",6, rotation=90))
annotate!(9.25, 1.25e4, text("PWR type SMR concepts",6, rotation=90))
annotate!(12.25, 1.25e4, text("HTR type SMR concepts",6, rotation=90))
annotate!(15.25, 1.25e4, text("SFR type SMR concepts",6, rotation=90))
savefig("fig-lcoe_mcs.pdf")