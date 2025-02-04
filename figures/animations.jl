using BIRD
using Makie
using CairoMakie
using LaTeXStrings

function wc_lengths()

    L = Observable(1.0)

    # Define morse potential
    lih = Morse(mA=1.0u"u", mB=6.941u"u", 
    ν=1405.0u"cm^-1", νχ=23.1679u"cm^-1", 
    re=1.595u"Å", De=2.641u"eV")
    
    # Most important transitions according to the sensitivity analysis
    overt_idx = [(21, 29), (22, 29), (23, 29), (25, 29)]
    overt_wvn = [get_transition_wvn(lih, l[1], l[2]) for l in overt_idx]

    fig = Figure(fontsize=18)
    ax = Axis(fig[1,1], yticks=0:0.5:1.5, ylabel=L"Relative Density of States $$", xlabel=L"Wavenumber (cm$^{-1}$)",
    title=@lift(title($L)))
    ax.xlabelsize=22
    ax.ylabelsize=22
    ylims!(ax, 0.0, 1.9)
    xlims!(ax, 300, 2020)
    ax.xticks = 500:500:2500
    hlines!(ax, [1.0], linestyle=:dash, color=:black)

    νrange = 300:5:6000

    ratio = @lift begin
        dos_at_L($L)
    end

    clrs = @lift begin
        get_colors($L, overt_wvn)
    end

    lines!(ax, νrange, ratio, linewidth=2)
    vlines!(ax, overt_wvn, linestyle=:dot, color=clrs, ymax=1.0, linewidth=3)

    # Add label to vertical lines
    for k in eachindex(overt_wvn)
        (i,j) = overt_idx[k]
        ν = overt_wvn[k] + 30
        text!(ax, ν, 0.06, text=L"%$i \rightarrow %$j", align=(:left, :top), fontsize=18, rotation=π/2)
    end

    Lmin = 2.0
    Lmax = 10.0
    step = 0.05
    Lvals = vcat(Lmin:step:Lmax, Lmax:-step:Lmin)
    record(fig, "dos.gif", Lvals;
        framerate = 30) do l
        L[] = l
    end
end

function dos_at_L(L)
    # Conversion factor
    μm2Å = 1e4

    νrange = 300:5:6000
    ratio = zeros(length(νrange))

    # Get relative DOS
    for i in eachindex(ratio)
        ω = ν2ω(νrange[i])
        ratio[i] = Dc(ω, L*μm2Å) / D0(ω)
    end

    return ratio
end

function get_colors(L, transitions)
    # Conversion factor
    μm2Å = 1e4
    return [(Dc(ν2ω(ν), L*μm2Å)/D0(ν2ω(ν)) > 1.0 ? :teal : :salmon3) for ν in transitions]
end

function title(L)
    lround = round(L, digits=1)
    L"Cavity Length = %$lround $\mu$m"
end