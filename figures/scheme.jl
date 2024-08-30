using BIRD
using Makie, WGLMakie, CairoMakie
using FileIO
using LaTeXStrings

function scheme()

    fig = Figure(size=(700, 400), fontsize=12)

    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1:2,1], xlabel=L"Bond Length (\AA)$$", ylabel=L"Energy (eV)$$", yticks=0:2:6)
    ax2 = Axis(gd[1,2], aspect = DataAspect())
    ax3 = Axis(gd[1,3], aspect=DataAspect())
    ax4 = Axis(gd[2,2], backgroundcolor=:grey34, xticks=1000:2000:9000, yticks=0.2:0.2:1.3, xlabel=L"$q$ (cm$^{-1}$)", ylabel=L"Energy (eV)$$")
    ax5 = Axis(gd[2,3], backgroundcolor=:grey34, xticks=1000:2000:9000, yticks=0.2:0.2:1.3, xlabel=L"$q$ (cm$^{-1}$)")

    # Hide decorations and spines
    #hidedecorations!(ax1)
    hidexdecorations!(ax1, label=false, ticks=false, ticklabels=false)
    hideydecorations!(ax1, label=false, ticks=false, ticklabels=false)
    hidespines!(ax1, :t, :r)
    hidedecorations!(ax2)
    hidespines!(ax2)
    hidedecorations!(ax3)
    hidespines!(ax3)
    hideydecorations!(ax5, ticks=false)

    # Plot morse potential
    plot_morse!(ax1)

    # Plot cavity figures
    img1 = FileIO.load(assetpath(joinpath(@__DIR__, "WCcavity.png")))
    img2 = FileIO.load(assetpath(joinpath(@__DIR__, "SCcavity.png")))
    image!(ax2, rotr90(img1))
    image!(ax3, rotr90(img2))

    # Plot dispersions
    dispersion_bare!(ax4)
    dispersion!(ax5)

    linkaxes!(ax4, ax5)
    ylims!(ax5, 0.1239, 0.9)
    xlims!(ax5, 0, 6000)

    colgap!(gd, 1, 0)
    colgap!(gd, 2, 20)

    fig

end

function plot_morse()
    fig = Figure()
    ax = Axis(fig[1,1])
    plot_morse!(ax)
    fig
end

function plot_morse!(ax)

    m = Morse(mA=1u"u", mB=19u"u", ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    # Morse energy
    E(r) = m.De * (1 - exp(-m.a*(r-m.re)))^2

    levels = [0, 1, 2, 23]
    ν2eV = 	1.23981e-4
    for l in vcat(0:1:22, 23)
        En = ν2eV*energy(m, l)
        rmax = m.re - (1/m.a) * log(1 - sqrt(En/m.De))
        rmin = m.re - (1/m.a) * log(1 + sqrt(En/m.De))
        lines!(ax, [rmin, rmax < 4 ? rmax : 4.0], [En, En], color=(:lightblue4, 0.5))
    end

    rvals = 0.6:0.01:4
    lines!(ax, rvals, [E(r) for r in rvals], color=:midnightblue, linewidth=4)

    # Arrows
    b1 = 11
    t1 = 15
    b2 = 9
    t2 = 5
    c = :salmon3
    x1 = 1.2
    xs = [x1, 0.85]
    ys = ν2eV .* ([energy(m, b1), energy(m, b2)])
    us = [0.0, 0.0]
    δ = 200
    vs = ν2eV .* ([get_transition_wvn(m, b1, t1)-δ, -get_transition_wvn(m, b2, t2)+1500])
    arrows!(ax, xs, ys, us, vs, linewidth = 3, color=c, arrowsize=15)
    text!(ax, x1+0.05, ν2eV*(energy(m,b1)+0.5*get_transition_wvn(m, b1, t1)), 
    text=L"k^\text{Abs}", align=(:left, :center), fontsize=22)

    text!(ax, 0.9, ν2eV*(energy(m,b2)-0.5*get_transition_wvn(m, b2, t2)), 
    text=L"k^\text{Em}", align=(:left, :center), fontsize=22)

    # kloss
    R = 1.0
    x0 = 1.8
    y0 = ν2eV*energy(m, 23)
    fk(x) = x > R ? y0 + R : y0 + √(R^2 - (x-R)^2)
    xvals = 0:0.01:2 

    #arrows!(ax, [1.8], [ν2eV*energy(m, 23)], [0.0], [0.8], linewidth=3, linestyle=:dot, color=c)
    lines!(ax, x0 .+ xvals, fk.(xvals), color=c, linewidth=3, linestyle=:dot)
    scatter!(ax, [x0 + xvals[end]], [fk(xvals[end])], marker=:rtriangle, color=c, markersize=16)
    text!(ax, 2.7, 6.7, text=L"k^\text{loss}", align=(:left, :center), fontsize=22)
end

function dispersion()
    fig = Figure()
    ax = Axis(fig[1,1], backgroundcolor=:grey17)
    dispersion!(ax)
    fig
end

function dispersion!(ax)

    # Conversion factor cm-1 to eV
    cm2eV = 0.0001239841984332003

    qvals = 10:100:7000
    L = 15.00e-4 # cm
    νM = 4000
    g = 700

    for m in 0:5
        q0 = m*π/L
        clp = [(:gold2, Pc_lp(q, m, L, νM, g)) for q in qvals]
        cup = [(:gold2, Pc_up(q, m, L, νM, g)) for q in qvals]
        lines!(ax, qvals, [cm2eV*UP_energy(q, m, L, νM, g) for q in qvals], linewidth=3, color=cup)
        lines!(ax, qvals, [cm2eV*LP_energy(q, m, L, νM, g) for q in qvals], linewidth=3, color=clp)
        lines!(ax, qvals, [cm2eV*sqrt(q^2 + q0^2) for q in qvals], linestyle=:dash, color=:white)
    end
    hlines!(ax, [νM*cm2eV], linestyle=:dash, color=:skyblue2)
    text!(ax, 2400, 0.25, text=L"\text{TM}_0", align=(:left, :top), color=:white)
end

function dispersion_bare!(ax)

    # Conversion factor cm-1 to eV
    cm2eV = 0.0001239841984332003
    qvals = 10:100:7000
    L = 15.00e-4 # cm
    νM = 4000*cm2eV

    for m in 0:5
        q0 = m*π/L
        lines!(ax, qvals, [cm2eV*sqrt(q^2 + q0^2) for q in qvals], linewidth=3, color=:gold2)
        lines!(ax, qvals, [cm2eV*sqrt(q^2 + q0^2) for q in qvals], linestyle=:dash, color=:white)
    end
    text!(ax, 2400, 0.25, text=L"\text{TM}_0", align=(:left, :top), color=:white)
    hlines!(ax, [νM], linestyle=:dash, color=:skyblue2)
end

function UP_energy(q, m, L, νM, g)

    νc = Ecav(q, m, L)
    ν2 = 0.5 * (νc^2 + νM^2 + g^2 + sqrt( (νM^2 + g^2 - νc^2)^2  + 4*νc^2 * g^2))

    return sqrt(ν2)
end

function LP_energy(q, m, L, νM, g)

    νc = Ecav(q, m, L)
    ν2 = 0.5 * (νc^2 + νM^2 + g^2 - sqrt( (νM^2 + g^2 - νc^2)^2  + 4*νc^2 * g^2))

    return sqrt(ν2)
end

function Ecav(q, m, L)
    q0 = m*π/L
    return sqrt(q^2 + q0^2)
end

function Pc_up(q, m, L, νM, g)
    νc = Ecav(q, m, L)
    ν = UP_energy(q, m, L, νM, g)
    Ω = (ν^2 - νM^2 - g^2)^2
    return Ω / (Ω + νc^2*g^2)
end

function Pc_lp(q, m, L, νM, g)
    νc = Ecav(q, m, L)
    ν = LP_energy(q, m, L, νM, g)
    Ω = (ν^2 - νM^2 - g^2)^2
    return Ω / (Ω + νc^2*g^2)
end