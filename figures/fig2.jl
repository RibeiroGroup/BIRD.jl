function plot_dos()

    # Create figure and axis
    fig = Figure(fontsize=18)
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j], yticks=_latexthis(0:0.5:1.5)) for i = 1:2, j = 1:2]
    linkaxes!(axs[1], axs[2])
    linkaxes!(axs[3], axs[4])

    ylims!(axs[2], -0.05, 1.9)
    xlims!(axs[2], 1000, 5000)
    axs[2].xticks=_latexthis(2000:1000:4000)

    ylims!(axs[4], -0.05, 1.9)
    xlims!(axs[4], 300, 2020)
    axs[4].xticks=_latexthis(500:500:2500)

    axs[1].title = "HF"
    axs[3].title = "LiH"

    # Hide y-axis bells and whistles in the right column
    for i in [3,4]
        hideydecorations!(axs[i], grid=false, ticks=false)
    end

    # Hide x-axis bells and whistles in the upper panels
    for i in [1, 3]
        hidexdecorations!(axs[i], grid=false, ticks=false)
    end

    # Adjust x labels
    for i in [2,4]
    end

    Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")
    Label(gd[1:2, 0], L"Relative Density of States $$", rotation=π/2)

    colgap!(gd, 1, 1.5)
    rowgap!(gd, 2, 1.5)

    ##################### LiH ##############################
    # Define morse potential
    lih = Morse(mA=1.0u"u", mB=6.941u"u", 
    ν=1405.0u"cm^-1", νχ=23.1679u"cm^-1", 
    re=1.595u"Å", De=2.641u"eV")

    path = joinpath(@__DIR__, "../QMcalcs/LiH/dip/dips.h5")
    r = vcat(0.0, h5read(path, "rvals")[4:end])
    dip = 0.529177 .* vcat(0.0, h5read(path, "dips")[4:end]) # Convert from e⋅bohr to e⋅Å
    # Get interpolation
    itp = interpolate((r,), dip, Gridded(Linear()))
    # Dipole function. Return 0 for out of bounds
    lih_dipole(x) = x > 6.6 ? 0.0 : itp(x)

    # Get matrix with transition energies
    E = transition_energy_matrix(lih)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(lih_dipole, lih)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = get_transport_matrix(E, V, DOS0, 2000)
    k0 = eigmin(J0)
    println("LiH")
    println(k0)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(21, 29), (22, 29), (23, 29), (25, 29)]
    overt_wvn = [get_transition_wvn(lih, l[1], l[2]) for l in overt_idx]

    wc_dos!(axs[3], 9.4, E, V, overt_wvn, overt_idx, k0=k0, T=2000, letter='c')
    wc_dos!(axs[4], 9.6, E, V, overt_wvn, overt_idx, k0=k0, T=2000, letter='d')

    ##################### HF ##############################
    # Define morse potential
    hf = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV")

    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4) 

    println(hf.nmax)
    # Get matrix with transition energies
    E = transition_energy_matrix(hf)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(hf_dipole, hf)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free-space rate
    J0 = get_transport_matrix(E, V, DOS0, 4000)
    k0 = eigmin(J0)
    println("HF")
    println(k0)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(16, 23), (15, 23), (14, 23), (17, 23), (19, 23), (13, 23)]
    overt_wvn = [get_transition_wvn(hf, l[1], l[2]) for l in overt_idx]

    wc_dos!(axs[1], 3.03, E, V, overt_wvn, overt_idx, k0=k0, T=4000, letter='a')
    wc_dos!(axs[2], 3.21, E, V, overt_wvn, overt_idx, k0=k0, T=4000, letter='b')

    fig
end

function wc_dos!(ax, L, E, V, transitions, transitions_labels; T=300, transition_label_y=0.02, k0, letter)

    # Conversion factor
    μm2Å = 1e4

    # Get cavity density of states (regime = 1)
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)

    # Get cavity rate
    Jc = get_transport_matrix(E, V, DOSc, T)
    kc = eigmin(Jc)
    println("L = $L kc = $kc")

    νrange = 300:5:12000
    ratio = zeros(length(νrange))

    # Plot relative DOS
    for i in eachindex(ratio)
        ω = ν2ω(νrange[i])
        ratio[i] = Dc(ω, L*μm2Å) / D0(ω)
    end

    clrs = [(Dc(ν2ω(ν), L*μm2Å)/D0(ν2ω(ν)) > 1.0 ? :teal : :salmon3) for ν in transitions]

    # Plot relative DOS
    lines!(ax, νrange, ratio, linewidth=2)

    # Vertical reference line
    hlines!(ax, [1.0], linestyle=:dash, color=:black)

    # Plot vertical lines for the desired transitions
    vlines!(ax, transitions, linestyle=:dot, color=clrs, ymax=0.85, linewidth=3)

    # Add label to vertical lines
    for k in eachindex(transitions)
        (i,j) = transitions_labels[k]
        ν = transitions[k] + 30
        text!(ax, ν, transition_label_y, text=L"%$i \rightarrow %$j", align=(:left, :top), fontsize=15, rotation=π/2)
    end

    kck0 = "$(round(kc/k0, digits=2))"
    text!(ax, 0.02, 0.93, text=L"(%$letter) $L_C =$ %$L $\mu$m $k_c/k_0 =$ %$kck0", align=(:left, :center), space=:relative, fontsize=15)
end