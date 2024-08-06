using Makie, WGLMakie, CairoMakie

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
    b = 8
    t = 14
    c = :tan1
    xs = [0.8, 1.3]
    ys = ν2eV .* ([energy(m, b), energy(m, t-2)])
    us = [0.0, 0.0]
    δ = 500
    vs = ν2eV .* ([get_transition_wvn(m, b, t)-δ, -get_transition_wvn(m, t+2, b)+1500])
    arrows!(ax, xs, ys, us, vs, linewidth = 3, color=c, arrowsize=15)
    arrows!(ax, [1.8], [ν2eV*energy(m, 23)], [0.0], [0.4], linewidth=3, Linestyle=:dash, color=c)
    text!(ax, 0.85, 4.45, text=L"k^\text{Abs}", align=(:left, :center), fontsize=28)
    text!(ax, 1.35, 3.9, text=L"k^\text{Em}", align=(:left, :center), fontsize=28)
    text!(ax, 1.85, 6.5, text=L"k^\text{loss}", align=(:left, :center), fontsize=28)
end

function dispersion()
    fig = Figure()
    ax = Axis(fig[1,1], backgroundcolor=:grey17)
    dispersion!(ax)
    fig
end

function dispersion!(ax)

    qvals = 10:100:7000
    L = 15.00e-4 # cm
    νM = 4000
    g = 700

    for m in 1:5
        q0 = m*π/L
        clp = [(:gold2, Pc_lp(q, m, L, νM, g)) for q in qvals]
        cup = [(:gold2, Pc_up(q, m, L, νM, g)) for q in qvals]
        lines!(ax, qvals, [UP_energy(q, m, L, νM, g) for q in qvals], linewidth=3, color=cup)
        lines!(ax, qvals, [LP_energy(q, m, L, νM, g) for q in qvals], linewidth=3, color=clp)
        lines!(ax, qvals, [sqrt(q^2 + q0^2) for q in qvals], linestyle=:dash, color=:white)
    end
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