using BIRD
using Makie
using WGLMakie
using LinearAlgebra
using LaTeXStrings

function find_halflife(J)

    t = 0.0
    N0 = zeros(size(J,1))
    N0[10] = 1.0

    δt = 0.001
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

    error("Couldn't find half-life in 100 steps of 0.1 s")

end

function perturbation_analysis(;δ=1.2, T=300)

    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
    m = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV", 
    μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    # Get matrix with transition energies
    E = transition_energy_matrix(m)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(m)

    # Get free space density of states (regime = 0)
    DOS = get_DOS_matrix(E, regime=0)

    K = find_pertbation_matrix(E, V, DOS, T=T, δ=δ)

    l, Kl = _linearize_offdiagonal(K)
    pe = reverse(sortperm(Kl))

    for k = 1:10
        (i,j) = l[pe][k]
        kenh = Kl[pe][k]
        println("$i -> $j  --- $kenh")
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

    # Get matrix with transition energies
    E = transition_energy_matrix(m)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(m)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, L=1.2,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          regime=1)

    J = get_transport_matrix(E, V, DOS0, 300)

    # Modify J such that the system becomes conservative
    J[1,end] = -J[end,end]

    pop0 = zeros(size(J,1))
    pop0[1] = 1.0
    pop1 = similar(pop0)

    δt = 0.01
    U = exp(-J*δt)

    Δ = 1.0
    iter = 1
    maxiter = 1000
    while Δ > 1e-10
        if iter > maxiter                                                                                                                                                                                                   
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


