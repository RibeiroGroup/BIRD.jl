function _latexthis(x)
    return (x, [L"%$s" for s in x])
end

function fig4(;T=300)

    axticks = [4000, 6000, 8000, 1000]
    ayticks = 1:2:15
    bxticks = [5600, 6000, 6400]
    byticks = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
    cxticks = [L"10^0", L"10^{-1}", L"10^{-2}", L"10^{-3}", L"10^{-4}"]
    ci = 15

    fig = Figure(resolution=(900, 400), fontsize=23)
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], xticks=_latexthis(axticks), yticks=_latexthis(ayticks), ylabel=L"Relative rate ($k_p/k_c$)")
    ax2 = Axis(gd[1,2], xticks=_latexthis(bxticks), yticks=_latexthis(byticks))
    ax3 = Axis(gd[1,3], xticks=(0:4, cxticks), yticks=_latexthis(ayticks))

    Label(gd[2,1:2], L"$\omega_M$ (cm$^{-1}$)")
    Label(gd[2,3], L"$\omega_M - \omega_{%$ci\rightarrow 23}$ (cm$^{-1}$)", tellwidth=false)

    xlims!(ax2, 5500, 6500)
    ylims!(ax2, 0.75, 1.3)

    ylims!(ax3, 1, 15.5)

    rowgap!(gd, 1, 8)
    colgap!(gd, 1, 10)
    colgap!(gd, 2, 10)

    # Add labels
    for (i, ax) in enumerate([ax1, ax2, ax3])
        text!(ax, 0.02, 0.99, text="($('`'+i))", align=(:left, :top), space=:relative, fontsize=20, font=:bold)
    end

    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
    m = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV", 
    μ0=7.27615208838288e-20u"C", ζ=0.081616u"Å^-4")

    # Highlighted overtone
    overt = get_transition_wvn(m, ci, 23)

    # Get matrix with transition energies
    E = transition_energy_matrix(m)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(m)

    # Get free space density of states (regime = 0)
    μm2Å = 1e4 # Conversion factor
    L = 1.08
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)

    # Get free space rate
    Jc = get_transport_matrix(E, V, DOSc, T)
    kc = eigmin(Jc)

    νM = 3000
    νMmax = 10000

    νMvals = Float64[]
    r = Float64[]
    while νM < νMmax

        # If νM is too close to a transition, add a small pertubation
        if minimum(abs.(E .- ν2ω(νM))) < 1e-3
            νM += 1e-2
            continue
        end

        DOSp = get_DOS_matrix(E, L=L*μm2Å, g=600, νM=νM, regime=2)
        Jp = get_transport_matrix(E, V, DOSp, T)
        kp = eigmin(Jp)

        push!(r, kp/kc)
        push!(νMvals, νM)

        νM += 0.5
    end

    hlines!(ax1, [1.0], color=:gray, linestyle=:dash)
    hlines!(ax2, [1.0], color=:gray, linestyle=:dash)

    lines!(ax1, νMvals, r)
    lines!(ax2, νMvals, r, linewidth=3)
    vlines!(ax2, [overt], color=:salmon3, linestyle=:dot, linewidth=2)

    # Saturation analysis
    prox = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

    gs = [100, 300, 600, 1000]

    for g in gs
        r = zeros(length(prox))
        for i in eachindex(prox)
            DOSp = get_DOS_matrix(E, L=L*μm2Å, g=g, νM=overt + 1/(10^(prox[i])), regime=2)
            Jp = get_transport_matrix(E, V, DOSp, T)
            kp = eigmin(Jp)
            r[i] = kp/kc
        end

        lines!(ax3, prox, r, label = L"%$g")
        scatter!(ax3, prox, r, label = L"%$g")
    end

    axislegend(L"$\Omega_R$ (cm$^{-1}$)", merge=true, position=:rb, labelsize=20)

    fig
end