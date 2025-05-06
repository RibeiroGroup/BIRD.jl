using BIRD
using Makie
using LinearAlgebra
using DataInterpolations
using Quadmath, DoubleFloats, ArnoldiMethod
using LaTeXStrings
using HDF5

function lowest_eigenval(M::Matrix{Float64})
    return eigmin(M)
end

function lowest_eigenval(M::Matrix{Float128})
    return partialschur(M, nev=1, which=:SR)[1].eigenvalues[1]
end

"""
    Quick validation for the number of bound states.
In this plot, we check that the predicted number of bound states (Morse)
matches what the numerov wavefunctions tell us. Bound wave functions must
vanish at r → ∞. 
"""
function plot_numerov_wvf()
    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    
    # Compute Numerov wavefunctions
    rgrid, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)

    plot_numerov_wvf(nali, rgrid, e, C, [53, 54, 55, 56])
end

function plot_numerov_wvf(M, r, e, C, nvals)

    De = potential(M, r[end])

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="R (Å)", ylabel="Energy (cm⁻¹)")

    println("Dissociation energy: $De cm⁻¹")

    for n in nvals
        println("n = $n | Energy: $(e[n+1])) cm⁻¹")
        lines!(ax, r, C[:,n+1], label="$n")
    end

    axislegend(ax, "n", position=:rt)

    fig
end

"""
    Plot kloss for each state 
"""
function plot_kloss()

    #### Get NaLi data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    
    # Compute Numerov wavefunctions
    rgrid, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)

    plot_kloss(nali, Dfunc, BIRD.D0, rgrid, e, C, 400)
end

function plot_kloss(M, dip_func, DOS_func, rgrid, e, C, T)

    # Free space decay rates
    kploss0 = BIRD.compute_kploss(dip_func, DOS_func, M.nmax, rgrid, e, C, T)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"Initial State ($n_i$)", ylabel=L"$k_\text{loss}$ (s$^{-1}$)")

    lines!(ax, 0:(length(kploss0)-1), kploss0)

    fig
end

function weak_coupling(;T=400, precision=Float128, Lvals=4:0.1:30)

    fig = Figure(fontsize=20, size=(700,600))
    gd = GridLayout(fig[1,1])
    ax1 = Axis(gd[1,1:2], xlabel=L"Cavity Length ($\mu$m)", ylabel=L"Relative rate ($k_c/k_0$)")
    ax2 = Axis(gd[2,1], ylabel=L"Relative Density of States $$")
    ax3 = Axis(gd[2,2], ylabel=L"Relative Density of States $$")

    Label(gd[3, 1:2], L"Wavenumber (cm$^{-1}$)")

    rowgap!(gd, 1, 5)
    rowgap!(gd, 2, 0.5)

    text!(ax1, 0.03, 0.93, text=L"(a)$$", align=(:left, :center), space=:relative, fontsize=20)

    xlims!(ax1, 4, 27)
    ylims!(ax1, 0.7, 1.7)

    #### Get NaLi data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    
    # Compute Numerov wavefunctions
    rgrid, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)

    # Free space decay rates
    kploss0 = BIRD.compute_kploss(Dfunc, BIRD.D0, nali.nmax, rgrid, e, C, T)

    # Get matrix with transition energies
    E = transition_energy_matrix(nali)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(Dfunc, nali, rvals=rgrid)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T, kploss=kploss0))
    k0 = lowest_eigenval(J0)

    kcr = zeros(length(Lvals))

    dos_Lvals = Float64[]
    dos_kvals = Float64[]

    println("Computing kc values...")
    for i in eachindex(Lvals)

        println("L = $(Lvals[i])")
        # Conversion factor
        μm2Å = 1e4

        # Get cavity density of states (regime = 1)
        DOSc = get_DOS_matrix(E, L=Lvals[i]*μm2Å, regime=1)

        # Update kploss (it also changes due to the the modified DOS!)
        kplossc = BIRD.compute_kploss(Dfunc, x->BIRD.Dc(x, Lvals[i]*μm2Å), nali.nmax, rgrid, e, C, T)

        # Get relative cavity rate
        J = precision.(get_transport_matrix(E, V, DOSc, T, kploss=kplossc))
        kcr[i] = lowest_eigenval(J) / k0

        if Lvals[i] == 10.5 || Lvals[i] == 15.5
            push!(dos_Lvals, Lvals[i])
            push!(dos_kvals, round(kcr[i], digits=2))
        end
    end
    println("Done!")

    # Plot horizontal line at 1.0
    hlines!(ax1, [1.0], color=:gray, linestyle=:dash)

    # Plot ratio line and scatter
    lines!(ax1, Lvals, kcr, color=:darkgreen)

    hidexdecorations!(ax2, ticks=false, ticklabels=false)
    hidexdecorations!(ax3, ticks=false, ticklabels=false)
    hideydecorations!(ax3, ticks=false)
    linkaxes!(ax2, ax3)
    xlims!(ax3, 80, 800)
    ylims!(ax3, 0.3, 2)

    println("Preparing DOS plots...")
    idx = [(12, 14), (13, 15), (11, 13), (14, 16), (15, 17), (10, 12), (16, 18), (9, 11), (17, 19), (18, 20)]
    transitions = [BIRD.ω2ν(E[i+1,j+1]) for (i,j) in idx]

    dos!(ax2, L=dos_Lvals[1], r=dos_kvals[1], transitions=transitions, letter='b')
    dos!(ax3, L=dos_Lvals[2], r=dos_kvals[2], transitions=transitions, letter='c')
    println(Lvals[14])
    println("Done.")

    fig
end

function dos!(ax; L, r=nothing, transitions=[], letter='a')

    # Range of wave
    νrange = 50:1:800

    # Conversion factor
    μm2Å = 1e4

    ratio = zeros(length(νrange))

    # Plot relative DOS
    for i in eachindex(νrange)
        ω = ν2ω(νrange[i])
        ratio[i] = Dc(ω, L*μm2Å) / D0(ω)
    end

    clrs = [(Dc(ν2ω(ν), L*μm2Å)/D0(ν2ω(ν)) > 1.0 ? :teal : :salmon3) for ν in transitions]

    # Plot relative DOS
    lines!(ax, νrange, ratio, linewidth=2)

    # Vertical reference line
    hlines!(ax, [1.0], linestyle=:dash, color=:black)

    # Plot vertical lines for the desired transitions
    vlines!(ax, transitions, linestyle=:dot, color=clrs, ymax=0.89, linewidth=2)

    text!(ax, 0.05, 0.93, text=L"(%$letter) $L_C =$ %$L $\mu$m $k_c/k_0 =$ %$r", align=(:left, :center), space=:relative, fontsize=20)
end

function top10(δ=1.5; number=10, kploss=true, printout=false, precision=Float64)

    # Get molecular data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    # Anharmonic constant was obtained from the dissociation energy
    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.3273u"cm^-1", 
        re=2.889u"Å", De=0.882388u"eV")

    E = transition_energy_matrix(nali)
    V = transition_dipole_matrix(Dfunc, nali, rvals=0.0:0.01:50, overtones=true)
    DOS = get_DOS_matrix(E, regime=0)    

    if kploss
        # Compute Numerov wavefunctions
        rgrid, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)

        # Free space decay rates
        kploss = BIRD.compute_kploss(Dfunc, BIRD.D0, nali.nmax, rgrid, e, C, 2000)
    end

    # Compute the sensitivity as S = k(after perturbation) / k(before perturbation) - 1
    S = find_pertbation_matrix(E, V, DOS, T=2000, δ=δ, kploss=kploss)

    N = size(S,1)
    r = 0:(N-1)
    labels = [(i,j) for i = r, j = r]
    m = [i < j for i = r, j = r]

    newS = S[m]
    idxs = labels[m]

    p = sortperm(newS, rev=true)

    outS = newS[p][1:number]
    outidx = idxs[p][1:number]

    if printout
        for k in 1:number
            i,j = outidx[k]
            s = outS[k]
            println("Trasition $i → $j: $s   / Energy = $(BIRD.ω2ν(E[i+1,j+1]))")
        end
    end

    return outidx, outS

end

function _linearize_offdiagonal(M)
    N = size(M,1)

    return labels[m], M[m]
end

function compare_J(L1, L2; T=400)

    #### Get NaLi data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    
    # Compute Numerov wavefunctions
    rgrid, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)

    # Get matrix with transition energies
    E = transition_energy_matrix(nali)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(Dfunc, nali, rvals=rgrid)

    μm2Å = 1e4
    ##### L1
    # Get cavity density of states (regime = 1)
    DOSc1 = get_DOS_matrix(E, L=L1*μm2Å, regime=1)

    kploss1 = BIRD.compute_kploss(Dfunc, x->BIRD.Dc(x, L1*μm2Å), nali.nmax, rgrid, e, C, T)
    J1 = get_transport_matrix(E, V, DOSc1, T, kploss=kploss1)
    # J1 = get_transport_matrix(E, V, DOSc1, T)
    # k1 = partialschur(Float128.(J1), nev=1, which=:SR)[1].eigenvalues[1] / k0

    ##### L2
    # Get cavity density of states (regime = 1)
    DOSc2 = get_DOS_matrix(E, L=L2*μm2Å, regime=1)

    kploss2 = BIRD.compute_kploss(Dfunc, x->BIRD.Dc(x, L2*μm2Å), nali.nmax, rgrid, e, C, T)
    J2 = get_transport_matrix(E, V, DOSc2, T, kploss=kploss2)
    # J2 = get_transport_matrix(E, V, DOSc2, T)


    fig =Figure()
    ax1 = Axis(fig[1,1])
    hm = heatmap!(ax1, J2 ./ J1)

    ax2 = Axis(fig[1,2])
    heatmap!(ax2, DOSc2 ./ DOSc1)

    ax3 = Axis(fig[2,1:2], xlabel=L"Morse state ($n$)", ylabel=L"Decay rate $k_\text{loss}$ (s$^{-1}$)")
    scatter!(ax3, 0:(length(kploss1)-1), kploss1, label = L"%$L1 \AA")
    scatter!(ax3, 0:(length(kploss2)-1), kploss2, label = L"%$L2 \AA")

    Colorbar(fig[:, end+1], hm)

    fig
end

function kploss_sensitivity(δ; T=400, precision=Float128, printout=false)

    #### Get NaLi data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    # Compute Numerov wavefunctions
    rgrid, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)
    # Free space decay rates
    kploss0 = BIRD.compute_kploss(Dfunc, BIRD.D0, nali.nmax, rgrid, e, C, T)

    # Get matrix with transition energies
    E = transition_energy_matrix(nali)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(Dfunc, nali, rvals=rgrid)
    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T, kploss=kploss0))
    k0 = nothing
    if eltype(J0) == Float64
        k0 = eigmin(J0)
    else
        k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    end
    println("k0 = $k0")

    Svec = zeros(length(kploss0))

    for i in eachindex(kploss0)
        kploss = zeros(length(kploss0))
        kploss .= kploss0
        kploss[i] = δ * kploss[i]

        J = precision.(get_transport_matrix(E, V, DOS0, T, kploss=kploss))
        if eltype(J) == Float64
            k = eigmin(J)
            Svec[i] = k / k0 - 1
        else
            k = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
            Svec[i] = k / k0 - 1
        end

        if printout
            if i == 1
                println("-----   Sensitivity Score for each decay to the continuum rate (δ = $δ) -----")
            end
            # println("Initial State: $i     Sensitivity: $(round(Svec[i], digits=5))")
            println("Initial State: $i     Sensitivity: $(Svec[i])")
        end
    end

    return Svec
end

function strong_coupling(fname="strong_coupling_ratios.h5"; mindetuning=1.0, labeln=3)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    r = h5read(fname, "ratios")
    νMvals = h5read(fname, "freqs")

    E = BIRD.ω2ν.(transition_energy_matrix(nali))

    m1 = r .!= 0
    m2 = [all(abs(x - e) > mindetuning for e in E) for x in νMvals]

    m = m1 .& m2

    ### Figure
    fig = Figure(size=(500,500), fontsize=20)
    ax = Axis(fig[1,1], ylabel=L"Relative rate ($k_c/k_0$)", xlabel=L"$\omega_M$ (cm$^{-1}$)")
    # ax2 = Axis(fig[2,1])
    # ax3 = Axis(fig[2,2])
    scatter!(ax, νMvals[m], r[m])
    lines!(ax, νMvals[m], r[m])

    # Reference lines for free-space and weak-coupling case
    hlines!(ax, [1.0], color=:grey20, linestyle=:dash)
    hlines!(ax, [1.5], color=:darkgreen, linestyle=:dash)

    # idx = [(12, 14), (13, 15), (11, 13), (14, 16), (15, 17), (10, 12), (16, 18), (9, 11), (17, 19), (18, 20)]
    # transitions = [get_transition_wvn(nali, i, j) for (i,j) in idx]

    # vlines!(transitions, linestyle=:dash, color=:black)

    # plot_sc_dos!(ax2, L=5, νM=300, gs=[25, 75, 125])
    # plot_sc_dos!(ax3, L=5, νM=400, gs=[25, 75, 125])


    fig
end


function generate_strong_coupling(outname; νMvals=200:0.5:500, L=5, T=400, precision=Float128, g=125)

    #### Get NaLi data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.33u"cm^-1", 
        re=2.895u"Å", De=0.88u"eV")

    
    # Compute Numerov wavefunctions
    rgrid, e, C = BIRD.get_continuum_states(nali, N=3000, rmin=0.0, rmax=50.0)

    # Free space decay rates
    kploss0 = BIRD.compute_kploss(Dfunc, BIRD.D0, nali.nmax, rgrid, e, C, T)

    # Get matrix with transition energies
    E = transition_energy_matrix(nali)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(Dfunc, nali, rvals=rgrid)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T, kploss=kploss0))
    k0 = lowest_eigenval(J0)
    μm2Å = 1e4

    νMvals = collect(νMvals)

    idx = [(12, 14), (13, 15), (11, 13), (14, 16), (15, 17), (10, 12), (16, 18), (9, 11), (17, 19), (18, 20), 
    (8, 10), (19, 21), (20, 22), (7, 9), (21, 23), (26, 29), (25, 28), (27, 30), (24, 27), (5, 6)]

    for (i,j) in idx
        push!(νMvals, get_transition_wvn(nali, i, j) + 0.1)
        push!(νMvals, get_transition_wvn(nali, i, j) - 0.1)
    end

    sort!(νMvals)

    r = zeros(length(νMvals))
    for (i,νM) in enumerate(νMvals)

        println("νM = $νM")

        # Get kploss under the polaritonic enviroment (νM)
        kploss = BIRD.compute_kploss(Dfunc, x->BIRD.Dp(x, L*μm2Å, ν2ω(νM), ν2ω(g)), nali.nmax, rgrid, e, C, T)

        # Compute density of states with polaritons characterized by νM
        DOSp = get_DOS_matrix(E, L=L*μm2Å, g=g, νM=νM, regime=2)
        if any(isinf, DOSp)
            println("Oops! That frequency blows up!")
            continue
        end

        # New transport matrix
        Jp = precision.(get_transport_matrix(E, V, DOSp, T, kploss=kploss))

        # Rate constant
        kp = lowest_eigenval(Jp)

        # Ratio (relative rate)
        r[i] = kp/k0
    end

    h5write(outname, "ratios", r)
    h5write(outname, "freqs", collect(νMvals))
end

function plot_sc_dos!(ax; L=26.8, νM=380, gs=[0, 75, 125])

    # Conversion factor
    μm2Å = 1e4

    # Plot strong coupling
    listyles = [:solid, :dash, :dot]
    n = 1
    for (g,l) in zip(gs, listyles)
        # Range o ν values. Increased density of points around νM
        stopgap = 4
        νrange = vcat(100:10:(νM-50),
        (νM-50):1:(νM-stopgap),
        (νM+stopgap):1(νM+50),
        (νM+50):10:2000)

        ratio = zeros(length(νrange))
        # Plot relative DOS
        for i in eachindex(ratio)
            ω = ν2ω(νrange[i])
            ratio[i] = Dp(ω, L*μm2Å, ν2ω(νM), ν2ω(g)) / D0(ω)
        end

        # Plot relative DOS
        lines!(ax, νrange, ratio, linewidth=3, label=L"%$g", linestyle=l)
        lines!(ax, νrange, ratio, linewidth=1, color=Makie.wong_colors()[n])
        n += 1
    end

    # Vertical reference line
    hlines!(ax, [1.0], linestyle=:dash, color=:black)
end