using BIRD
using Makie
using WGLMakie
using LinearAlgebra
using LaTeXStrings

function wc_k_ratio_vs_L(;T=300, Lvals = 0.5:0.01:5, circle=true, Lcircles = [1.05, 1.08, 3.16, 3.31])

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

    # Get free-space rate
    J0 = get_transport_matrix(E, V, DOS0, T)
    k0 = eigmin(J0)
    println(k0)

    kcr = zeros(length(Lvals))

    for i in eachindex(Lvals)
        # Conversion factor
        μm2Å = 1e4

        # Get cavity density of states (regime = 1)
        DOSc = get_DOS_matrix(E, L=Lvals[i]*μm2Å, regime=1)

        # Get relative cavity rate
        J = get_transport_matrix(E, V, DOSc, T,)
        kcr[i] = eigmin(J) / k0
    end

    # Create figure and axis
    fig = Figure(fontsize=25, resolution=(600,400))
    ax1 = Axis(fig[1,1], xlabel=L"Cavity Length ($\mu$m)", ylabel=L"Relative rate ($k_c/k_0$)", xticks=0:1:Lvals[end], xminorticks=0.5:1:4.5, xminorticksvisible = true)

    # Plot horizontal line at 1.0
    hlines!(ax1, [1.0], color=:gray, linestyle=:dash)

    # Plot ratio line and scatter
    lines!(ax1, Lvals, kcr)
    scatter!(ax1, Lvals, kcr, markersize=5)


    # Adjust limits of the plot
    xmin, xmax = Lvals[1], Lvals[end]
    ymin, ymax = 0.9, 1.35
    ylims!(ax1, ymin, ymax)
    xlims!(ax1, xmin, xmax)


    # Create new axis to draw highlights on top 
    axcirc = Axis(fig[1,1])
    hidedecorations!(axcirc)
    xlims!(axcirc, 0,1)
    ylims!(axcirc, 0, 1)

    # Draw circles
    if circle
        for L in Lcircles
            i = findfirst(x->x==L, Lvals)
            k = kcr[i]

            # Find the coordinates in the ax1 system and convert it to axcirc
            Y = (k - ymin) / (ymax - ymin)
            X = (L - xmin)  / (xmax - xmin)
            arc!(axcirc, Point2f(X, Y), 0.03, -π, π, color=:red, linewidth=4)
        end
    end

    fig
end
