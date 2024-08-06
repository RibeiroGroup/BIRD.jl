using Plots
using BIRD
using LinearAlgebra
using LaTeXStrings

function plot_dipole(;SI=false)

    rvals = 0:0.01:4

    m = Morse(mA=1u"u", mB=19u"u", ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", re=0.926u"Å", De=6.121646u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    dips = zeros(length(rvals))

    for i in eachindex(dips)
        dips[i] = BIRD.get_dipole(m, rvals[i])
    end

    if SI
        e = 1.60217663e-19
        Å2m = 1e-10
        plot(rvals, 1e30 * dips * e * Å2m, label="")
        ylabel!("μ (C⋅m 10^-30)")
    else
        plot(rvals, dips, label="")
        ylabel!("μ (e⋅Å)")
    end

    xlabel!("R (Å)")
end

function plot_transition_dipole(;SI=false)

    vvals = 0:23

    m = Morse(mA=1u"u", mB=19u"u", ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", re=0.926u"Å", De=6.121646u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    dips = zeros(length(vvals))

    for i in eachindex(dips)
        v = vvals[i]
        dips[i] = BIRD.transition_dipole(m, v,v)
    end

    if SI
        e = 1.60217663e-19
        Å2m = 1e-10
        plot(vvals, 1e30 * dips * e * Å2m, label="")
        scatter!(vvals, 1e30 * dips * e * Å2m, label="")
        ylabel!("⟨i|μ|j⟩ (C⋅m 10^-30)")
    else
        plot(vvals, dips, label="")
        scatter!(vvals, dips, label="")
        ylabel!("μ (e⋅Å)")
        ylabel!("⟨i|μ|j⟩ (e⋅Å)")
    end

    xlabel!("v")
end

function plot_morse(v)

    rvals = 0:0.01:8
    m = Morse(mA=1u"u", mB=19u"u", ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", re=0.926u"Å", De=6.121646u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")
    E = BIRD.energy(m, v) / 8065.61
    wfn = BIRD.compute_wfn(m, v, rvals)

    V = [m.De*(1 - exp(-m.a*(r - m.re)))^2 for r in rvals]

    plot(rvals, V, label="Potential")
    plot!(rvals, wfn .+ E, label="Wavefunction")
    ylims!(-1,7)

end

function average_r()
end

function testrate()
    N = zeros(24)
    N[1] = 1.0

    m = Morse(mA=1u"u", mB=19u"u", ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", re=0.926u"Å", De=6.121646u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")
    J = BIRD.get_transport_matrix(m, 300, kloss=1e10)
    tvals = 0:0.01:2
    out = zeros(length(tvals))

    for i in eachindex(tvals)
        t = tvals[i]
        U = exp(-J*t)
        new = U*N
        out[i] = sum(new)
    end

    plot(tvals, out)
end

function ratio()
    m = Morse(mA=1u"u", mB=19u"u", ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", re=0.926u"Å", De=6.121646u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")
    Lvals = 0.5:1:80 # μm -> must convert to Å before using it
    r = zeros(length(Lvals))

    for i in eachindex(Lvals)
        L = Lvals[i]
        r[i] = eigmin(BIRD.get_cav_transport_matrix(m, L*10000, 300)) / eigmin(BIRD.get_transport_matrix(m, 300))
        println(r[i])
    end

    plot(Lvals, r)
    xlabel!("L (μm)")
    ylabel!("k_d/k_d^0")
end

function dispersion()
    qvals = 0:100:50000

    L = 1 # μm
    mvals = 1:3

    plot()
    for m in mvals
        ωcs = [BIRD.cav_wvn(q, L/10000, m) for q in qvals]
        plot!(qvals, ωcs, label="Cavity (m=$m)", linestyle=:dot)

        UP = [BIRD.UP_freq(ωc, 6000, 100) for ωc in ωcs]
        LP = [BIRD.LP_freq(ωc, 6000, 100) for ωc in ωcs]

        plot!(qvals, UP, label="UP$m")
        plot!(qvals, LP, label="LP$m")
    end

    plot!()
end

function plot_Ds(;νM=6050, L=10000, g=200)

    nuvals = 1000:100:12000

    D0s = BIRD.D0.(BIRD.ν2ω.(nuvals))

    plot(nuvals, [BIRD.Dc(BIRD.ν2ω(nu), 10000) for nu in nuvals] ./ D0s, label="Weak")
    plot!(nuvals, [BIRD.Dp(BIRD.ν2ω(nu), 10000, BIRD.ν2ω(νM), BIRD.ν2ω(g)) for nu in nuvals] ./ D0s, label="Strong")

    xlabel!("Wavenumber (cm^-3)")
    ylabel!("Relative Density of States (D/D0)")
    plot!()
end