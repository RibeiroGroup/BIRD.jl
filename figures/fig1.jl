using LinearAlgebra
using Interpolations
using HDF5

# Turn ticks into latex ticks using LaTeXStrings
function _latexticks(x)
    return (x, [L"%$s" for s in x])
end

function wc_k_ratio_vs_L(; circle=true)

    ### HF
    # Create figure and axis
    fig = Figure(fontsize=25, size=(600,400))
    ax1 = Axis(fig[1,1], xlabel=L"Cavity Length ($\mu$m)", ylabel=L"Relative rate ($k_c/k_0$)", 
    xticks=_latexticks(0:2:20), xminorticks=0.5:1:4.5, xminorticksvisible = true,
    yticks=_latexticks(0.9:0.1:1.3))

    # Adjust limits of the plot
    Lvals = 0.5:0.01:20
    xmin, xmax = Lvals[1], Lvals[end]
    ymin, ymax = 0.8, 1.35
    ylims!(ax1, ymin, ymax)
    xlims!(ax1, xmin, xmax)

    # Define a Morse potential. Parameters are taken from Kaluza and Muckerman 1993 (https://doi.org/10.1063/1.466305)
    # Note that μ0 from the reference has been converted from e.s.u to C using 1 e.s.u = 3.335640951982e-10 C
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
    println("HF Free-space dissociation rate: ", k0)

    kcr = zeros(length(Lvals))
    for i in eachindex(Lvals)
        # Conversion factor
        μm2Å = 1e4

        # Get cavity density of states (regime = 1)
        DOSc = get_DOS_matrix(E, L=Lvals[i]*μm2Å, regime=1)

        # Get relative cavity rate
        J = get_transport_matrix(E, V, DOSc, 4000)
        kcr[i] = eigmin(J) / k0
    end

    # Plot horizontal line at 1.0
    hlines!(ax1, [1.0], color=:gray, linestyle=:dash)

    # Plot ratio line and scatter
    lines!(ax1, Lvals, kcr, label=L"HF$$")

    # Create new axis to draw highlights on top 
    axcirc = Axis(fig[1,1])
    hidedecorations!(axcirc)
    xlims!(axcirc, 0, 1)
    ylims!(axcirc, 0, 1)

    # Draw circles 
    Lcircles = [3.03, 3.31]
    if circle
        for L in Lcircles
            i = findfirst(x->x==L, Lvals)
            k = kcr[i]

            # Find the coordinates in the ax1 system and convert it to axcirc
            Y = (k - ymin) / (ymax - ymin)
            X = (L - xmin)  / (xmax - xmin)
            arc!(axcirc, Point2f(X, Y), 0.01, -π, π, color=:red, linewidth=2)
        end
    end

    ############################################# LiH
    lih = Morse(mA=1.0u"u", mB=6.941u"u", 
    ν=1405.0u"cm^-1", νχ=23.1679u"cm^-1", 
    re=1.595u"Å", De=2.641u"eV")
    println(lih.nmax)

    # Read data for dipole function
    path = joinpath(@__DIR__, "../QMcalcs/LiH/dip/dips.h5")
    r = vcat(0.0, h5read(path, "rvals"))
    dip = 0.529177 .* vcat(0.0, h5read(path, "dips")) # Convert from e⋅bohr to e⋅Å
    # Get interpolation
    itp = interpolate((r,), dip, Gridded(Linear()))
    # Dipole function. Return 0 for out of bounds
    lih_dipole(x) = x > 6.6 ? 0.0 : itp(x)

    # Get matrix with transition energies
    E = transition_energy_matrix(lih)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(lih_dipole, lih, rvals=0:0.01:25)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free-space rate
    J0 = get_transport_matrix(E, V, DOS0, 2000)
    k0 = eigmin(J0)
    println("LiH Free-space dissociation rate: ", k0)

    kcr = zeros(length(Lvals))
    for i in eachindex(Lvals)
        # Conversion factor
        μm2Å = 1e4

        # Get cavity density of states (regime = 1)
        DOSc = get_DOS_matrix(E, L=Lvals[i]*μm2Å, regime=1)

        # Get relative cavity rate
        J = get_transport_matrix(E, V, DOSc, 2000)
        kcr[i] = eigmin(J) / k0
    end

    # Plot ratio line and scatter
    lines!(ax1, Lvals, kcr, label=L"LiH$$")

    # Draw circles 
    Lcircles = [9.4, 9.6]
    if circle
        for L in Lcircles
            i = findfirst(x->x==L, Lvals)
            k = kcr[i]

            # Find the coordinates in the ax1 system and convert it to axcirc
            Y = (k - ymin) / (ymax - ymin)
            X = (L - xmin)  / (xmax - xmin)
            arc!(axcirc, Point2f(X, Y), 0.01, -π, π, color=:red, linewidth=2)
        end
    end

    axislegend(ax1)
    fig
end