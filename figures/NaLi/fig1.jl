using BIRD
using DataInterpolations
using LinearAlgebra
using Quadmath, DoubleFloats, ArnoldiMethod

# Produces an upper panel with the Weak coupling Relative rate versus cavity length (L)
# Lower panels Relative density of states (weak coupling) versus wavenumber at two different values of L
function fig1(T; precision=Float64)

    # Create figure and axis
    fig = Figure(fontsize=20, size=(700,600))
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1:2], xlabel=L"Cavity Length ($\mu$m)", ylabel=L"Relative rate ($k_c/k_0$)")
    axs = [Axis(gd[2,1], ylabel=L"Relative Density of States $$"), Axis(gd[2,2])]

    La, Lb = 26.3, 26.8
    #La, Lb = 15.7, 26.6
    #La, Lb = 16.0, 37.8 
    #La, Lb = 18.8, 22.8 

    Lvals= 15.0:0.1:40
    ymin, ymax = 0.95, 1.30
    xmin, xmax = Lvals[1], Lvals[end]
    ylims!(ax1, ymin, ymax)
    xlims!(ax1, xmin, xmax)

    Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")

    # Create new axis to draw highlights on top 
    axcirc = Axis(gd[1,1:2])
    hidedecorations!(axcirc)
    xlims!(axcirc, 0, 1)
    ylims!(axcirc, 0, 1)

    linkaxes!(axs...)
    ylims!(axs[2], -0.3, 2.5)

    # Hide y-axis bells and whistles in the right column
    hideydecorations!(axs[2], grid=false, ticks=false)

    # Fit Morse
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    E = transition_energy_matrix(nali)
    V = transition_dipole_matrix(Dfunc, nali, rvals=0.0:0.01:50, overtones=true)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T))
    if eltype(J0) == Float64
        k0 = eigmin(J0)
    else
        k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    end
    println("NaLi ($T K)")
    println("k0 = ", k0)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(42, 54), (43, 54), (48, 54)]
    overt_wvn = [get_transition_wvn(nali, l[1], l[2]) for l in overt_idx]

    # Plot density of states in the lower panels
    wc_dos!(axs[1], La, E, V, 50:1:600, overt_wvn, overt_idx, k0=k0, T=T, letter='b', precision=precision)
    wc_dos!(axs[2], Lb, E, V, 50:1:600, overt_wvn, overt_idx, k0=k0, T=T, letter='c', precision=precision)

    # Plot ratios in the upper panel
    wc_k_ratios!(ax1, axcirc, E, V, T, Lvals, [La, Lb], xmin, xmax, ymin, ymax, "", circle=false, precision=precision)
    text!(ax1, 0.02, 0.93, text=L"$$(a)", align=(:left, :center), space=:relative, fontsize=20)

    rowgap!(gd, 1, 5)
    rowgap!(gd, 2, 0.5)
    fig
end

function wc_k_ratios!(ax, axcirc, E, V, T, Lvals, Lcircles, xmin, xmax, ymin, ymax, label; circle=true, precision=Float64)
    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free-space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T))
    if eltype(J0) == Float64 
        k0 = eigmin(J0)
    else
        k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    end
    println("$label Free-space dissociation rate: ", k0)

    kcr = zeros(length(Lvals))
    for i in eachindex(Lvals)
        # Conversion factor
        μm2Å = 1e4

        # Get cavity density of states (regime = 1)
        DOSc = get_DOS_matrix(E, L=Lvals[i]*μm2Å, regime=1)

        # Get relative cavity rate
        J = precision.(get_transport_matrix(E, V, DOSc, T))
        if eltype(J) == Float64
            kcr[i] = eigmin(J) / k0
        else
            kcr[i] = partialschur(J, nev=1, which=:SR)[1].eigenvalues[1] / k0
        end
    end

    # Plot horizontal line at 1.0
    hlines!(ax, [1.0], color=:gray, linestyle=:dash)

    # Plot ratio line and scatter
    lines!(ax, Lvals, kcr, label=label, color=:darkgreen)

    # Draw circles 
    if circle
        for L in Lcircles[1:2]
            i = findfirst(x->x==L, Lvals)
            k = kcr[i]

            # Find the coordinates in the ax1 system and convert it to axcirc
            Y = (k - ymin) / (ymax - ymin)
            X = (L - xmin)  / (xmax - xmin)
            arc!(axcirc, Point2f(X, Y), 0.01, -π, π, color=:red, linewidth=2)
        end
    end
end

function wc_dos!(ax, L, E, V, νrange, transitions, transitions_labels; T=300, transition_label_y=-0.2, k0, letter, precision=Float64)

    # Conversion factor
    μm2Å = 1e4

    # Get cavity density of states (regime = 1)
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)

    # Get cavity rate
    Jc = precision.(get_transport_matrix(E, V, DOSc, T))
    if eltype(Jc) == Float64
        kc = eigmin(Jc)
    else
        kc = partialschur(Jc, nev=1, which=:SR)[1].eigenvalues[1]
    end
    println("L = $L kc = $kc")

    xlims!(ax, νrange[1], νrange[end])
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
        ν = transitions[k] + 10
        text!(ax, ν, transition_label_y, text=L"%$i \rightarrow %$j", align=(:left, :top), fontsize=15, rotation=π/2)
    end

    kck0 = "$(round(Float64(kc/k0), digits=2))"
    text!(ax, 0.02, 0.93, text=L"(%$letter) $L_C =$ %$L $\mu$m $k_c/k_0 =$ %$kck0", align=(:left, :center), space=:relative, fontsize=18)
end

function LiH_fig2(T; precision=Float64)

    println("Start")
    # Create figure and axis
    fig = Figure(fontsize=18, size=(600,600))
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1:2], xlabel=L"Cavity Length ($\mu$m)", ylabel=L"Relative rate ($k_c/k_0$)")
    axs = [Axis(gd[2,1], ylabel=L"Relative Density of States $$", xticks=500:500:2000), Axis(gd[2,2], xticks=500:500:2000)]

    #La, Lb = 37.8, 38.30 
    La, Lb = 6.8, 6.9

    Lvals= 2.0:0.1:15
    ymin, ymax = 0.8, 1.5
    xmin, xmax = Lvals[1], Lvals[end]
    ylims!(ax1, ymin, ymax)
    xlims!(ax1, xmin, xmax)

    Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")

    # Create new axis to draw highlights on top 
    axcirc = Axis(gd[1,1:2])
    hidedecorations!(axcirc)
    xlims!(axcirc, 0, 1)
    ylims!(axcirc, 0, 1)

    linkaxes!(axs...)
    ylims!(axs[2], -0.3, 2.5)

    # Hide y-axis bells and whistles in the right column
    hideydecorations!(axs[2], grid=false, ticks=false)

    pes = get_pes("LiH")

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(pes.E, regime=0)
    println("A")

    # Get free space rate
    J0 = precision.(get_transport_matrix(pes.E, pes.V, DOS0, T))
    if eltype(J0) == Float64
        k0 = eigmin(J0)
    else
        k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    end
    println("LiH ($T K)")
    println(k0)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(18, 22), (17, 22), (20, 22)]
    overt_wvn = [get_transition_wvn(pes, l[1], l[2]) for l in overt_idx]

    wc_dos!(axs[1], La, pes.E, pes.V, 400:10:2200, overt_wvn, overt_idx, k0=k0, T=T, letter='a', precision=precision)
    wc_dos!(axs[2], Lb, pes.E, pes.V, 400:10:2200, overt_wvn, overt_idx, k0=k0, T=T, letter='b', precision=precision)

    wc_k_ratios!(ax1, axcirc, pes.E, pes.V, T, Lvals, [La, Lb], xmin, xmax, ymin, ymax, "", circle=true, precision=precision)

    fig
end

function NaH_fig2(T; precision=Float64)

    println("Start")
    # Create figure and axis
    fig = Figure(fontsize=18, size=(600,600))
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1:2], xlabel=L"Cavity Length ($\mu$m)", ylabel=L"Relative rate ($k_c/k_0$)")
    axs = [Axis(gd[2,1], ylabel=L"Relative Density of States $$", xticks=500:500:2000), Axis(gd[2,2], xticks=500:500:2000)]

    #La, Lb = 37.8, 38.30 
    La, Lb = 8.4, 8.5

    Lvals= 2.0:0.1:15
    ymin, ymax = 0.8, 1.5
    xmin, xmax = Lvals[1], Lvals[end]
    ylims!(ax1, ymin, ymax)
    xlims!(ax1, xmin, xmax)

    Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")

    # Create new axis to draw highlights on top 
    axcirc = Axis(gd[1,1:2])
    hidedecorations!(axcirc)
    xlims!(axcirc, 0, 1)
    ylims!(axcirc, 0, 1)

    linkaxes!(axs...)
    ylims!(axs[2], -0.3, 2.5)

    # Hide y-axis bells and whistles in the right column
    hideydecorations!(axs[2], grid=false, ticks=false)

    pes = get_pes("NaH")

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(pes.E, regime=0)
    println("A")

    # Get free space rate
    J0 = precision.(get_transport_matrix(pes.E, pes.V, DOS0, T))
    if eltype(J0) == Float64
        k0 = eigmin(J0)
    else
        k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    end
    println("LiH ($T K)")
    println(k0)

    # Most important transitions according to the sensitivity analysis
    overt_idx = [(17, 21), (16, 22), (19, 22)]
    overt_wvn = [get_transition_wvn(pes, l[1], l[2]) for l in overt_idx]

    wc_dos!(axs[1], La, pes.E, pes.V, 400:10:2200, overt_wvn, overt_idx, k0=k0, T=T, letter='a', precision=precision)
    wc_dos!(axs[2], Lb, pes.E, pes.V, 400:10:2200, overt_wvn, overt_idx, k0=k0, T=T, letter='b', precision=precision)

    wc_k_ratios!(ax1, axcirc, pes.E, pes.V, T, Lvals, [La, Lb], xmin, xmax, ymin, ymax, "", circle=true, precision=precision)
    fig
end