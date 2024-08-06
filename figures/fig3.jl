using BIRD
using Makie
using WGLMakie
using LinearAlgebra
using LaTeXStrings

function plot_sc_dos(;L=1.08, νM=6000, g=200, T=300)

    # Conversion factor
    μm2Å = 1e4

    # Adjust x labels
    xtcks = 0:2000:15000
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
    xlims!(ax, 2000, 11000)

    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
    m = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV", 
    μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    # Plot strong coupling
    for g in [0, 200, 600, 1200]
        # Range o ν values. Increased density of points around νM
        stopgap = 4
        νrange = vcat(500:10:(νM-50),
        (νM-50):1:(νM-stopgap),
        (νM+stopgap):1(νM+50),
        (νM+50):10:12000)

        ratio = zeros(length(νrange))
        # Plot relative DOS
        for i in eachindex(ratio)
            ω = ν2ω(νrange[i])
            ratio[i] = Dp(ω, L*μm2Å, ν2ω(νM), ν2ω(g)) / D0(ω)
        end

        # Plot relative DOS
        lines!(ax, νrange, ratio, linewidth=3, label=L"%$g")
    end

    # Vertical reference line
    hlines!(ax, [1.0], linestyle=:dash, color=:black)
    axislegend(L"$\Omega_R$ (cm$^{-1}$)")

    fig
end