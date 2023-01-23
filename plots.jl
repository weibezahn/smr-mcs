##### plots #####
# requires results from the main file

using CairoMakie

##### boxplots Monte Carlo simulation results #####
# requires results for all 15 reactor concepts

function mcs_plot(mcs_results, title::String, ylabel::String)

    xticks_wc = names(mcs_results)[1:9];
    xticks_ht = names(mcs_results)[10:12];
    xticks_sf = names(mcs_results)[13:15];

    mcs_boxplot = Figure();

    ax_wc = Axis(mcs_boxplot[1,1], xticks = (1:length(xticks_wc), xticks_wc), ylabel = ylabel);
    ax_wc.xticklabelrotation = π / 3;
    ax_wc.yticklabelrotation = π / 2;
    ax_wc.xticklabelalign = (:right, :center);

    for i in 1:length(xticks_wc)
        boxplot!(ax_wc, fill(i,n), mcs_results[!,i], color = :green)
    end;

    Label(mcs_boxplot[1, 1, Top()], "BWR & PWR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    ax_ht = Axis(mcs_boxplot[1,2], xticks = (1:length(xticks_ht), xticks_ht));
    ax_ht.xticklabelrotation = π / 3;
    ax_ht.yticklabelrotation = π / 2;
    ax_ht.xticklabelalign = (:right, :center);

    for i in 1:length(xticks_ht)
        boxplot!(ax_ht, fill(i,n), mcs_results[!,i+length(xticks_wc)], color = :red)
    end;

    Label(mcs_boxplot[1, 2, Top()], "HTR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    ax_sf = Axis(mcs_boxplot[1,3], xticks = (1:length(xticks_sf), xticks_sf));
    ax_sf.xticklabelrotation = π / 3;
    ax_sf.yticklabelrotation = π / 2;
    ax_sf.xticklabelalign = (:right, :center);

    for i in 1:length(xticks_sf)
        boxplot!(ax_sf, fill(i,n), mcs_results[!,i+length(xticks_wc)+length(xticks_sf)], color = :orange)
    end;

    Label(mcs_boxplot[1, 3, Top()], "SFR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    Label(mcs_boxplot[0, :], title, fontsize = 24, font = "Noto Sans Bold", color = (:black, 0.25))

    colsize!(mcs_boxplot.layout, 1, Relative(6/10));
    colsize!(mcs_boxplot.layout, 2, Relative(2/10));
    colsize!(mcs_boxplot.layout, 3, Relative(2/10));

    return mcs_boxplot

end

fig_mcs_npv = mcs_plot(npv_results, "NPV", "[USD/MW]")
fig_mcs_lcoe = mcs_plot(lcoe_results, "LCOE", "[USD/MWh]")

save("fig-mcs_npv.pdf", fig_mcs_npv)
save("fig-mcs_lcoe.pdf", fig_mcs_lcoe)

##### heatmaps sensitivity indices #####

function si_plot(si_results, title::String)

    si_s = filter(:si => index -> index == "S", si_results);
    si_st = filter(:si => index -> index == "ST", si_results);
    xticks = si_s.var;
    yticks = names(si_results)[3:end];
    data_s = Matrix(si_s[:,3:end]);
    data_st = Matrix(si_st[:,3:end]);

    si_heatmap = Figure();

    ax_s = Axis(si_heatmap[1, 1], xticks = (1:length(xticks), xticks), yticks = (1:length(yticks), yticks));
    ax_s.xticklabelrotation = π / 3;
    ax_s.xticklabelalign = (:right, :center);
    hmap_s = heatmap!(ax_s, data_s, colormap = :deep, colorrange = (0, 1));
    for i in 1:length(xticks), j in 1:length(yticks)
        txtcolor = data_s[i, j] > 0.15 ? :white : :black
        text!(ax_s, "$(round(data_s[i,j], digits = 2))", position = (i, j),
            color = txtcolor, fontsize = 12, align = (:center, :center))
    end
    Label(si_heatmap[1, 1, Top()], "first-oder effect", font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    Colorbar(si_heatmap[1, 2], hmap_s; width = 15, ticksize = 15);

    ax_st = Axis(si_heatmap[1, 3], xticks = (1:length(xticks), xticks), yticks = (1:length(yticks), yticks), yaxisposition = :right);
    ax_st.xticklabelrotation = π / 3;
    ax_st.xticklabelalign = (:right, :center);
    hmap_st = heatmap!(ax_st, data_st, colormap = :deep, colorrange = (0, 1));
    for i in 1:length(xticks), j in 1:length(yticks)
        txtcolor = data_st[i, j] > 0.15 ? :white : :black
        text!(ax_st, "$(round(data_st[i,j], digits = 2))", position = (i, j),
            color = txtcolor, fontsize = 12, align = (:center, :center))
    end
    Label(si_heatmap[1, 3, Top()], "total-oder effect", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    Label(si_heatmap[0, :], title, fontsize = 24, font = "Noto Sans Bold", color = (:black, 0.25));

    return si_heatmap

end

fig_si_npv = si_plot(si_npv_results, "NPV")
fig_si_lcoe = si_plot(si_lcoe_results, "LCOE")

save("fig-si_npv.pdf", fig_si_npv)
save("fig-si_lcoe.pdf", fig_si_lcoe)