using BIRD
using Makie
using WGLMakie
using LinearAlgebra
using LaTeXStrings

# Turn ticks into latex ticks using LaTeXStrings
function _latexthis(x)
    return (x, [L"%$s" for s in x])
end

function find_halflife(m, T)

    J = get_transport_matrix(m, T)

    t = 0.0
    N0 = [boltzmann_fac(m, 0, i, T) for i in 0:m.nmax]
    N0 = N0 ./ sum(N0)
    println(N0)

    δt = 1e6
    U = exp(-δt*J)

    N = U*N0
    for i in 2:1000
        N = U*N
        P = sum(N)
        #println("t = $(i*δt)  P = $P")
        if P < 0.5
            return i*δt
        end
    end

    error("Couldn't find half-life in 100 steps of $δt s")

end

function perturbation_analysis(; δ=1.2, T=300)
    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
    m = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV")

    perturbation_analysis(m; δ=δ, T=300)
end

function perturbation_analysis(dip, m; δ=1.2, T=300)

    # Get matrix with transition energies
    E = transition_energy_matrix(m)

    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4) 
    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(dip, m)

    # Get free space density of states (regime = 0)
    DOS = get_DOS_matrix(E, regime=0)

    K = find_pertbation_matrix(E, V, DOS, T=T, δ=δ)

    l, Kl = _linearize_offdiagonal(K)
    pe = reverse(sortperm(Kl))

    for k = 1:20
        (i,j) = l[pe][k]
        kenh = Kl[pe][k]
        wvn = BIRD.get_transition_wvn(m, i, j)
        println("$i -> $j  --- $kenh  || $wvn")
    end

    fig, ax, hm = heatmap(0:23, 0:23, K .- 1)
    Colorbar(fig[:, end+1], hm)
    fig
end

function hmap(M)
    fig, ax, hm = heatmap(0:23, 0:23, M)
    Colorbar(fig[:, end+1], hm)
    fig
end

function find_pertbation_matrix(E, V, DOS; δ=1.2, T=300)

    # Get free space rate
    J0 = get_transport_matrix(E, V, DOS, T)
    k0 = eigmin(J0)
    #k0 = log(2)/(find_halflife(J0))

    N = size(E,1)

    K = similar(E)
    K .= 0.0

    for i in 1:N
        K[i,i] = 1.0
        for j in (i+1):N
            DOS[i,j] *= δ
            DOS[j,i] *= δ
            J = get_transport_matrix(E, V, DOS, T)
            K[i,j] = eigmin(J) / k0
            K[j,i] = K[i,j]
            DOS[i,j] /= δ
            DOS[j,i] /= δ
        end
    end

    return K
end

function _linearize_offdiagonal(M)
    N = size(M,1)
    r = 0:(N-1)
    labels = [(i,j) for i = r, j = r]
    m = [i < j for i = r, j = r]

    return labels[m], M[m]
end

function _linearize_offdiagonal(M1, M2)
    N = size(M,1)
    r = 0:(N-1)
    labels = [(i,j) for i = r, j = r]
    m = [i < j for i = r, j = r]

    return labels[m], M1[m], M2[m]
end

function find_eq_composition()

    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
    m = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV", 
    μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    find_eq_composition(m)
end

function find_eq_composition(m)

    # Get matrix with transition energies
    E = transition_energy_matrix(m)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(m)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    J = get_transport_matrix(E, V, DOS0, 300, kloss=0.0)

    # Modify J such that the system becomes conservative
    #J[1,end] = -J[end,end]

    #pop0 = zeros(size(J,1))
    #pop0[1] = 1.0
    pop0 = [boltzmann_fac(m, 0, i, 300) for i in 0:m.nmax]
    pop0 = pop0 ./ sum(pop0)
    println(pop0)
    pop1 = similar(pop0)

    δt = 1
    U = exp(-J*δt)

    Δ = 1.0
    iter = 1
    maxiter = 10000
    while Δ > 1e-10
        if iter > maxiter                                                                                                                                                                                                   
            println(pop1)
            error("Couldnt find a steady state population in $maxiter iterations with δt = $δt")
        end

        pop1 = U*pop0

        # Renormalize
        pop1 = pop1 ./ sum(pop1)

        Δ = sum(abs.(pop1 .- pop0))

        pop0 = pop1
        iter += 1
    end

    return J, pop1
end

function path_diagram()

    # Calculations
    J, pop = find_eq_composition()

    rates = similar(J)
    for i = axes(rates,1)
        for j = axes(rates, 2)
            rates[i,j] = J[i,j]*pop[j]
        end
    end

    l,r = _linearize_offdiagonal(rates)

    r = abs.(r) ./ maximum(abs.(r))

    # Figure
    fig = Figure(resolution=(400,1000))
    ax = Axis(fig[1,1])
    hidexdecorations!(ax)
    hideydecorations!(ax, label=false, ticks=false, ticklabels=false)


    hlines!(ax, collect(0:23), xmin=0.2, xmax=0.8, color=:black, linewidth=4)
    xlims!(ax, 0, 10)

    # Plot fundamental transitions
    #for i = eachindex(pop)
    #    x = 5
    #    y = i-1

    #    u = 0
    #    v = 0.8

    #    t = i != 24 ? fundrates[i] : 1.0

    #    arrows!(ax, [x], [y], [u], [v], linewidth=3, color=(:salmon3, t))
    #end
    for k = eachindex(r)
        (i,j) = l[k]
        if i == 0
            println(j)
        end
        if i == 0 && j == 23
            continue
        end
        x = 2 + (j-i)*0.5
        y = i + 0.2
        u = 0
        v = (j-i) - 0.4

        c = i == j - 1 ? :teal : :salmon3

        arrows!(ax, [x], [y], [u], [v], linewidth=3, color=(c, r[k]))
    end
    #for i in eachindex(pop)
    #    if i == 24
    #        continue
    #    end
    #    Jrow = pop[i] .* J[:,i]

    #    _, j = findmax(abs, Jrow)

    #    x = 2 + i*0.2
    #    y = i-1
    #    u = 0
    #    v = j - i + 1
    #    arrows!(ax, [x], [y], [u], [v], linewidth=3, color=(:salmon3, 1.0))
    #end

    fig

end

function plot_wvf(m, n)

    fig = Figure()
    ax = Axis(fig[1,1])

    rvals = 0.0:0.01:20
    lines!(ax, rvals, [potential(m, r) for r in rvals], color=:midnightblue, linewidth=4)

    ylims!(ax, -100, 10000) 

    psi = 5000 .* compute_wfn(m, n, rvals) .+ energy(m, n)
    

    lines!(ax, rvals, psi)
    fig

end

function boltzmann_fac(m, n1, n2, T)

    E1 = 2π * BIRD.ħ * BIRD.c_cm * energy(m, n1)
    E2 = 2π * BIRD.ħ * BIRD.c_cm * energy(m, n2)

    kT = T * BIRD.k

    exp(-E2/kT) / exp(-E1/kT)
end

function planks_dist(T)
    ħ = BIRD.ħ
    k = BIRD.k
    nuvals = 100:50:10000
    out = zeros(length(nuvals))

    for i in eachindex(out)
        ω = BIRD.ν2ω(nuvals[i])
        out[i] = D0(ω) * ħ * ω * (1 / (exp(ħ*ω/(k*T)) - 1))
    end

    lines(nuvals, out)
end

function plot_kd(nuvals=[1000.0, 1500, 2000], Des=[0.2, 0.5, 0.8], T=1000)

    fig = Figure()
    ax = Axis(fig[1,1])
    for De in Des
        ks = zeros(length(nuvals))

        i = 1
        for ν in nuvals
            anh = ν^2 / (4* 8065.54393734921 * De)
            m = BIRD.Morse(mA=40u"u", mB=20.0u"u", ν=ν*1u"cm^-1", νχ= anh*1u"cm^-1", re=0.926u"Å",
            De=De*1u"eV", μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

            E = transition_energy_matrix(m)
            DOS = get_DOS_matrix(E)
            V = transition_dipole_matrix(m, rvals=0:0.01:50)
            J = get_transport_matrix(E, V, DOS, T)

            ks[i] = eigmin(J)
            println(m.nmax)
            i += 1
        end

        lines!(ax, nuvals, ks)
        scatter!(ax, nuvals, ks)
    end

    fig
end

function compare_dipoles()


    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4) 

    # Read data for dipole function
    path = joinpath(@__DIR__, "../QMcalcs/MgNe/dip/ZCxyz/dips.h5")
    r = vcat(0.0, h5read(path, "rvals"))
    println(r)
    dip = 0.529177 .* vcat(0.0, h5read(path, "dips")) # Convert from e⋅bohr to e⋅Å
    itp1 = interpolate((r,), dip, Gridded(Linear()))
    mgne_dipole(x) = x > 19 ? 0.0 : itp1(x)


    path = joinpath(@__DIR__, "../QMcalcs/LiH/dip/dips.h5")
    r = vcat(0.0, h5read(path, "rvals")[6:end])
    println(r)
    dip = 0.529177 .* vcat(0.0, h5read(path, "dips")[6:end]) # Convert from e⋅bohr to e⋅Å
    itp2 = interpolate((r,), dip, Gridded(Linear()))
    lih_dipole(x) = x > 6.6 ? 0.0 : itp2(x)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="R (Å)", ylabel="Dipole (e⋅Å)")

    rvals = (0.0:0.01:10)
    lines!(ax, rvals, hf_dipole.(rvals), label="HF")
    lines!(ax, rvals, lih_dipole.(rvals), label="LiH")
    lines!(ax, rvals, mgne_dipole.(rvals), label="MgNe2+")

    axislegend(labelsize=20, titlesize=20, position=:rt)
    fig
end