function plot_sc_dos(;L=26.8, νM=380, gs=[0, 75, 125])

    # Conversion factor
    μm2Å = 1e4

    # Adjust x labels
    xtcks = 0:100:1000
    xlticks= [L"%$x" for x in xtcks]

    ytcks = 0:0.5:2.0
    ylticks= [L"%$y" for y in ytcks]

    # Create figure and axis
    fig = Figure(fontsize=25)
    gd = fig[1,1] = GridLayout()
    ax = Axis(gd[1,1], xticks=(xtcks, xlticks), yticks=(ytcks, ylticks), 
    ylabel=L"Relative Density of States $$",
    xlabel=L"Wavenumber (cm$^{-1}$)")

    ylims!(ax, -0.05, 2.1)
    xlims!(ax, 100, 810)

    # Plot strong coupling
    listyles = [:solid, :dash, :dot]
    n = 1
    for (g,l) in zip(gs, listyles)
        # Range o ν values. Increased density of points around νM
        stopgap = 4
        νrange = vcat(100:10:(νM-50),
        (νM-50):1:(νM-stopgap),
        (νM+stopgap):1(νM+50),
        (νM+50):10:2000)

        ratio = zeros(length(νrange))
        # Plot relative DOS
        for i in eachindex(ratio)
            ω = ν2ω(νrange[i])
            ratio[i] = Dp(ω, L*μm2Å, ν2ω(νM), ν2ω(g)) / D0(ω)
        end

        # Plot relative DOS
        lines!(ax, νrange, ratio, linewidth=3, label=L"%$g", linestyle=l)
        lines!(ax, νrange, ratio, linewidth=1, color=Makie.wong_colors()[n])
        n += 1
    end

    # Vertical reference line
    hlines!(ax, [1.0], linestyle=:dash, color=:black)
    axislegend(L"$\Omega_R$ (cm$^{-1}$)", labelsize=20, titlesize=20, position=:rt)

    fig
end