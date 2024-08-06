using BIRD
using Makie
using WGLMakie
using LinearAlgebra
using LaTeXStrings

function plot_dos(L=1.2; T=300)

    # Create figure and axis
    fig = Figure(fontsize=23)
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j], xticks=0:2000:15000) for i = 1:2, j = 1:2]
    linkaxes!(axs...)

    ylims!(axs[end], 0, 1.9)
    xlims!(axs[end], 2000, 11000)

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
        tcks = 0:2000:15000
        axs[i].xticks=(tcks, string.(tcks))
    end

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
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = get_transport_matrix(E, V, DOS0, T)
    k0 = eigmin(J0)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(15, 23), (16, 23), (14, 23), (17, 23), (13, 23)]
    overt_wvn = [get_transition_wvn(m, l[1], l[2]) for l in overt_idx]

    Lvals = [1.05, 3.16, 1.08, 3.31]
    #Lvals = [0.54, 0.65, 0.8, 0.9]
    letters = ['a', 'c', 'b', 'd']
    for i in eachindex(Lvals)
        L = Lvals[i]
        ax = axs[i]
        letter = letters[i]
        wc_dos!(ax, L, E, V, overt_wvn, overt_idx)
        text!(ax, 0.02, 0.93, text=L"(%$letter) $L_C =$ %$L $\mu$m", align=(:left, :center), space=:relative)
    end

    #Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")
    Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")
    Label(gd[1:2, 0], L"Relative Density of States $$", rotation=π/2)

    colgap!(gd, 1, 1.5)
    rowgap!(gd, 2, 1.5)

    fig
end

function wc_dos!(ax, L, E, V, transitions, transitions_labels; T=300, transition_label_y=0.25)


    # Conversion factor
    μm2Å = 1e4

    # Get cavity density of states (regime = 1)
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)

    # Get cavity rate
    Jc = get_transport_matrix(E, V, DOSc, T)
    kc = eigmin(Jc)
    println("L = $L kc = $kc")

    νrange = 500:10:12000
    ratio = zeros(length(νrange))

    # Plot relative DOS
    for i in eachindex(ratio)
        ω = ν2ω(νrange[i])
        ratio[i] = Dc(ω, L*μm2Å) / D0(ω)
    end

    clrs = [(Dc(ν2ω(ν), L*μm2Å)/D0(ν2ω(ν)) > 1.0 ? :teal : :salmon3) for ν in transitions]

    # Plot relative DOS
    lines!(ax, νrange, ratio, linewidth=3)

    # Vertical reference line
    hlines!(ax, [1.0], linestyle=:dash, color=:black)

    # Plot vertical lines for the desired transitions
    vlines!(ax, transitions, linestyle=:dot, color=clrs, ymax=0.91, linewidth=4)

    # Add label to vertical lines
    for k in eachindex(transitions)
        (i,j) = transitions_labels[k]
        ν = transitions[k]
        text!(ax, ν, transition_label_y, text=L"%$i \rightarrow %$j", align=(:center, :bottom), fontsize=15, rotation=π/2)
    end
end