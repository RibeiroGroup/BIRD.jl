function fig4()
    fig = Figure(size=(800, 800), fontsize=20)
    gd = fig[1,1] = GridLayout()

    hf_axticks = [3500, 4000, 4500]
    hf_ayticks = 0:10:60
    hf_cxticks = [L"10^0", L"10^{-1}", L"10^{-2}", L"10^{-3}", L"10^{-4}"]
    hf_ci = 17

    # Fist axis
    hf_ax1 = Axis(gd[1,1], xticks=_latexthis(hf_axticks), yticks=_latexthis(hf_ayticks), ylabel=L"Relative rate ($k_p/k_0$)",
    xlabel=L"$\omega_M$ (cm$^{-1}$)")

    # Second axis
    hf_ax2 = Axis(gd[1,2], xticks=(0:4, hf_cxticks), yticks=_latexthis(hf_ayticks), xlabel=L"$\omega_M - \omega_{%$hf_ci\rightarrow 23}$ (cm$^{-1}$)")
    hideydecorations!(hf_ax2, grid=false, ticks=false)

    # Inset axis
    hf_inax = Axis(gd[1,1], width=Relative(0.65), height=Relative(0.65), halign=0.8, valign=0.8, backgroundcolor=:lightgray,
    xticks=_latexthis([3400, 3500, 3600]), yticks=_latexthis([0.8, 1.0, 1.2, 1.4, 1.6]),
    xticklabelsize=15, yticklabelsize=15)

    sc_hf!(hf_ax1, hf_ax2, hf_inax, hf_ci)

    lih_axticks = 1000:100:4000
    lih_ayticks = 0:5:30
    lih_cxticks = [L"10^0", L"10^{-1}", L"10^{-2}", L"10^{-3}", L"10^{-4}"]
    lih_ci = 22

    # Fist axis
    lih_ax1 = Axis(gd[2,1], xticks=_latexthis(lih_axticks), yticks=_latexthis(lih_ayticks), ylabel=L"Relative rate ($k_p/k_0$)",
    xlabel=L"$\omega_M$ (cm$^{-1}$)")

    # Second axis
    lih_ax2 = Axis(gd[2,2], xticks=(0:4, lih_cxticks), yticks=_latexthis(lih_ayticks), xlabel=L"$\omega_M - \omega_{%$lih_ci\rightarrow 29}$ (cm$^{-1}$)")
    hideydecorations!(lih_ax2, grid=false, ticks=false)

    # Inset axis
    lih_inax = Axis(gd[2,1], width=Relative(0.65), height=Relative(0.65), halign=0.8, valign=0.8, backgroundcolor=:lightgray,
    xticks=_latexthis([1350, 1400, 1450]), yticks=_latexthis([0.8, 1.0, 1.2, 1.4, 1.6]),
    xticklabelsize=15, yticklabelsize=15)

    sc_lih!(lih_ax1, lih_ax2, lih_inax, lih_ci)

    # Add labels
    for (i, ax) in enumerate([hf_ax1, hf_ax2, lih_ax1, lih_ax2])
        text!(ax, 0.02, 0.99, text="($('`'+i))", align=(:left, :top), space=:relative, fontsize=20, font=:bold)
    end

    fig
end

function sc_hf!(ax1, ax2, inax, ci;T=4000)

    # Adjust axes
    xlims!(ax1, 3350, 4600)
    ylims!(ax1, 0, 55)
    ylims!(ax2, 0, 55)

    insxmin = 3400
    insxmax = 3600
    insymin = 0.8
    insymax = 1.8
    xlims!(inax, insxmin, insxmax)
    ylims!(inax, insymin, insymax)

    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
    hf = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV")

    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4) 

    # Highlighted overtone
    overt = get_transition_wvn(hf, ci, 23)

    # Get matrix with transition energies
    E = transition_energy_matrix(hf)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(hf_dipole, hf)

    # Get free space density of states (regime = 0)
    μm2Å = 1e4 # Conversion factor
    L = 3.31
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = get_transport_matrix(E, V, DOS0, T)
    k0 = eigmin(J0)
    println("k0: ", k0)

    # Get cavity rate
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)
    Jc = get_transport_matrix(E, V, DOSc, T)
    kc = eigmin(Jc)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(15, 23), (16, 23), (14, 23), (17, 23), (13, 23)]
    overt_wvn = [get_transition_wvn(hf, l[1], l[2]) for l in overt_idx] .+ 1e-4

    νMvals = vcat(
        3300:1:3499,
        overt_wvn[4],
        3500:1.87:4685,
        overt_wvn[2],
        4695:5:6050,
        overt_wvn[1],
        6055:5:7590,
        overt_wvn[3],
        7595:7.3:9300,
        overt_wvn[5],
        9305:5:10000
    )

    r = zeros(length(νMvals))

    for (i,νM) in enumerate(νMvals)

        DOSp = get_DOS_matrix(E, L=L*μm2Å, g=400, νM=νM, regime=2)
        if any(isinf, DOSp)
            println(νM)
        end
        Jp = get_transport_matrix(E, V, DOSp, T)
        kp = eigmin(Jp)

        r[i] = kp/k0
    end

    hlines!(ax1, [1.0], color=:gray, linestyle=:dash)
    hlines!(inax, [1.0], color=:gray, linestyle=:dash)
    hlines!(ax1, [kc/k0], color=:goldenrod, linestyle=:dash)
    hlines!(inax, [kc/k0], color=:goldenrod, linestyle=:dash)

    lines!(ax1, νMvals, r)
    lines!(inax, νMvals, r)

    # Saturation analysis
    prox = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

    gs = [100, 200, 400]

    n = 1
    ms = [:circle, :rect, :utriangle]
    for g in gs
        r = zeros(length(prox))
        for i in eachindex(prox)
            DOSp = get_DOS_matrix(E, L=L*μm2Å, g=g, νM=overt + 1/(10^(prox[i])), regime=2)
            Jp = get_transport_matrix(E, V, DOSp, T)
            kp = eigmin(Jp)
            r[i] = kp/k0
        end

        lines!(ax2, prox, r, label = L"%$g")
        scatter!(ax2, prox, r, label = L"%$g", marker=ms[n])
        n += 1
    end

    vlines!(inax, [overt], linestyle=:dot, color=:firebrick)
    text!(inax, overt, 0.85, text=L"%$ci \rightarrow 23", align=(:left, :top), fontsize=15, rotation=π/2)

    axislegend(ax2, L"$\Omega_R$ (cm$^{-1}$)", merge=true, position=:rb, labelsize=20)
end

function sc_lih!(ax1, ax2, inax, ci;T=2000)

    # Adjust axes
    xlims!(ax1, 1350, 1770)
    ylims!(ax1, 0, 25)
    ylims!(ax2, 0, 25)
    insxmin = 1350
    insxmax = 1475
    insymin = 0.7
    insymax = 1.8
    xlims!(inax, insxmin, insxmax)
    ylims!(inax, insymin, insymax)

    # Define morse potential
    lih = Morse(mA=1.0u"u", mB=6.941u"u", 
    ν=1405.0u"cm^-1", νχ=23.1679u"cm^-1", 
    re=1.595u"Å", De=2.641u"eV")

    path = joinpath(@__DIR__, "../QMcalcs/LiH/dip/dips.h5")
    r = vcat(0.0, h5read(path, "rvals"))
    dip = 0.529177 .* vcat(0.0, h5read(path, "dips")) # Convert from e⋅bohr to e⋅Å
    # Get interpolation
    itp = interpolate((r,), dip, Gridded(Linear()))
    # Dipole function. Return 0 for out of bounds
    lih_dipole(x) = x > 6.6 ? 0.0 : itp(x)

    # Highlighted overtone
    overt = get_transition_wvn(lih, ci, 29)
    println(overt)

    # Get matrix with transition energies
    E = transition_energy_matrix(lih)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(lih_dipole, lih)

    # Get free space density of states (regime = 0)
    μm2Å = 1e4 # Conversion factor
    L = 9.6
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = get_transport_matrix(E, V, DOS0, T)
    k0 = eigmin(J0)
    println(k0)

    # Get cavity rate
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)
    Jc = get_transport_matrix(E, V, DOSc, T)
    kc = eigmin(Jc)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(21, 29), (22, 29), (23, 29), (25, 29)]
    overt_wvn = [get_transition_wvn(lih, l[1], l[2]) for l in overt_idx]

    νMvals = vcat(
        500:1:4000
    )

    r = zeros(length(νMvals))

    for (i,νM) in enumerate(νMvals)

        DOSp = get_DOS_matrix(E, L=L*μm2Å, g=200, νM=νM, regime=2)
        if any(isinf, DOSp)
            println(νM)
        end
        Jp = get_transport_matrix(E, V, DOSp, T)
        kp = eigmin(Jp)

        r[i] = kp/k0
    end

    hlines!(ax1, [1.0], color=:gray, linestyle=:dash)
    hlines!(inax, [1.0], color=:gray, linestyle=:dash)
    hlines!(ax1, [kc/k0], color=:goldenrod, linestyle=:dash)
    hlines!(inax, [kc/k0], color=:goldenrod, linestyle=:dash)

    lines!(ax1, νMvals, r)
    lines!(inax, νMvals, r)

    # Saturation analysis
    prox = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

    gs = [50, 100, 200]

    n = 1
    ms = [:circle, :rect, :utriangle]
    for g in gs
        r = zeros(length(prox))
        for i in eachindex(prox)
            DOSp = get_DOS_matrix(E, L=L*μm2Å, g=g, νM=overt + 1/(10^(prox[i])), regime=2)
            Jp = get_transport_matrix(E, V, DOSp, T)
            kp = eigmin(Jp)
            r[i] = kp/k0
        end

        lines!(ax2, prox, r, label = L"%$g")
        scatter!(ax2, prox, r, label = L"%$g", marker=ms[n])
        n += 1
    end

    vlines!(inax, [overt], linestyle=:dot, color=:firebrick)
    text!(inax, overt, 0.75, text=L"%$ci \rightarrow 29", align=(:left, :top), fontsize=15, rotation=π/2)

    axislegend(ax2, L"$\Omega_R$ (cm$^{-1}$)", merge=true, position=:rb, labelsize=20)
end