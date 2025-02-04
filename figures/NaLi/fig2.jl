using BIRD
using DataInterpolations
using LinearAlgebra
using Quadmath, DoubleFloats, ArnoldiMethod


function fig2(;T=2000, precision=Float64)

    #### Get NaLi data
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

    ### Initialize Figure
    fig = Figure(size=(900, 500), fontsize=20)
    gd = fig[1,1] = GridLayout()

    ci = 42
    #### Panel (a) and inset
    ax1 = Axis(gd[1,1], xlabel=L"$\omega_M$ (cm$^{-1}$)", ylabel=L"Relative rate ($k_p/k_0$)", xticks=300:100:900)
    inax = Axis(gd[1,1], width=Relative(0.55), height=Relative(0.65), halign=0.8, valign=0.8, backgroundcolor=:grey90,
    xticklabelsize=15, yticklabelsize=15)
    kp_vs_ωM!(ax1=ax1, inax=inax, mol=nali, ci=ci, E=E, V=V, T=T, L=26.8, g=200, precision=precision)

    #### Panel (b)
    ax2 = Axis(gd[1,2], xticks=1:2:9, xlabel=L"$\omega_M - \omega_{%$ci\rightarrow 54}$ (cm$^{-1}$)")
    saturation_analysis!(ax=ax2, mol=nali, i=ci, E=E, V=V, T=T, L=26.8, precision=precision, gs=[100, 150, 200])

    # Link and adjust axis
    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2, ticks=false, grid=false)

    xlims!(ax1, 300, 1000)
    xlims!(inax, 355, 410)
    ylims!(inax, 0.8, 2.5)
    ylims!(ax2, 0.5, 46)

    colsize!(gd, 1, Relative(2/3))

    # Add labels
    text!(ax1, 0.02, 0.99, text="(a)", align=(:left, :top), space=:relative, fontsize=20, font=:bold)
    text!(ax2, 0.02, 0.99, text="(b)", align=(:left, :top), space=:relative, fontsize=20, font=:bold)

    fig
end

function saturation_analysis(;T=2000, precision=Float64)
    #### Get NaLi data
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

    ### Initialize Figure
    fig = Figure(size=(900, 500), fontsize=20)
    gd = fig[1,1] = GridLayout()

    #### Panel (b)
    ci = 42
    ax2 = Axis(gd[1,1], xticks=1:2:9, xlabel=L"$\omega_M - \omega_{%$ci\rightarrow 54}$ (cm$^{-1}$)")
    saturation_analysis!(ax=ax2, mol=nali, i=ci, E=E, V=V, T=T, L=26.8, gs=[100, 150, 200],  precision=precision)

    fig
end

function saturation_analysis!(;ax, mol, i, E, V, T, L, gs=[75, 100, 125], precision=Float64)

    # Transition index i -> j
    j = mol.nmax
    ν = get_transition_wvn(mol, i, j)

    μm2Å = 1e4 # Conversion factor
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T))
    k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    # k0 = eigmin(J0)
    println("k0 = $(real(k0))")

    prox = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    n = 1
    ms = [:circle, :rect, :utriangle]
    maxy = 0.0
    for g in gs
        r = zeros(length(prox))
        for i in eachindex(prox)
            DOSp = get_DOS_matrix(E, L=L*μm2Å, g=g, νM = ν + prox[i], regime=2)
            Jp = precision.(get_transport_matrix(E, V, DOSp, T))
            kp = partialschur(Jp, nev=1, which=:SR)[1].eigenvalues[1]
            #kp = eigmin(Jp)
            r[i] = kp/k0
        end

        maxy = max(maxy, maximum(r))
        lines!(ax, prox, r, label = L"%$g")
        scatter!(ax, prox, r, label = L"%$g", marker=ms[n])
        n += 1
    end
    println("Maximum polariton enhancement for $i -> $j ($ν cm⁻¹) at T = $T: $maxy")

    axislegend(ax, L"$\Omega_R$ (cm$^{-1}$)", merge=true, position=:rt, labelsize=20)
end

function kp_vs_ωM!(;ax1, inax, mol, ci, E, V, T, L, g=125, precision=Float64)

    # Highlighted overtone
    overt = get_transition_wvn(mol, ci, 54)

    # Get free space density of states (regime = 0)
    μm2Å = 1e4 # Conversion factor
    DOS0 = get_DOS_matrix(E, regime=0)

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS0, T))
    k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    println("k0 = $k0")

    # Get weak coupling cavity rate
    DOSc = get_DOS_matrix(E, L=L*μm2Å, regime=1)
    Jc = precision.(get_transport_matrix(E, V, DOSc, T))
    kc = partialschur(Jc, nev=1, which=:SR)[1].eigenvalues[1]
    println("kc = $kc")

    νMvals = sort(vcat(
        300:10:350, 355:1:378, 386:1:409, 410:10:1000,
        [get_transition_wvn(mol, 39, 54)+i for i in vcat(-5:0.1:-1, 1:0.1:5)],
        [get_transition_wvn(mol, 40, 54)+i for i in vcat(-5:0.1:-1, 1:0.1:5)],
        [get_transition_wvn(mol, 41, 54)+i for i in vcat(-5:0.1:-1, 1:0.1:5)],
        [get_transition_wvn(mol, 42, 54)+i for i in vcat(-5:0.1:-1, 1:0.1:5)],
        [get_transition_wvn(mol, 43, 54)+i for i in vcat(-5:0.1:-1, 1:0.1:5)],
    ))

    r = zeros(length(νMvals))
    for (i,νM) in enumerate(νMvals)

        DOSp = get_DOS_matrix(E, L=L*μm2Å, g=g, νM=νM, regime=2)
        if any(isinf, DOSp)
            println(νM)
        end
        Jp = precision.(get_transport_matrix(E, V, DOSp, T))
        kp = partialschur(Jp, nev=1, which=:SR)[1].eigenvalues[1]

        r[i] = kp/k0
    end

    hlines!(ax1, [1.0], color=:grey20, linestyle=:dash)
    hlines!(inax, [1.0], color=:grey20, linestyle=:dash, linewidth=2)
    hlines!(ax1, [kc/k0], color=:darkgreen, linestyle=:dash)
    hlines!(inax, [kc/k0], color=:darkgreen, linestyle=:dash, linewidth=2)

    lines!(ax1, νMvals, r)
    lines!(inax, νMvals, r)
    # scatter!(inax, νMvals, r)

    vlines!(inax, [overt], linestyle=:dot, color=:firebrick)
    text!(inax, overt, 1.25, text=L"%$ci \rightarrow 54", align=(:left, :top), fontsize=17, rotation=π/2)
end

function maximum_pol_effects(Tvals;precision=Float64)

    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.3273u"cm^-1", 
        re=2.889u"Å", De=0.882388u"eV")

    E = transition_energy_matrix(nali)
    V = transition_dipole_matrix(Dfunc, nali, rvals=0.0:0.01:50, overtones=true)

    # kploss
    #rvals, e, C = BIRD.get_continuum_states(nali, N=2000, rmin=0.0, rmax=50.0)

    fig = Figure(fontsize=20)
    ax = Axis(fig[1,1], xlabel=L"$i \rightarrow 54\;$ overtone resonant to $\omega_M$", ylabel=L"k_p/k_0", xticks=0:5:54, title="Maximum BIRD enhancement")

    for T in Tvals
        DOS0 = get_DOS_matrix(E, regime=0)
        #kploss0 = BIRD.compute_kploss(Dfunc, BIRD.D0, nali.nmax, rvals, e, C, 1000)
        #J0 = precision.(get_transport_matrix(E, V, DOS0, T, kploss=kploss0))
        J0 = precision.(get_transport_matrix(E, V, DOS0, T))
        k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]

        x = zeros(nali.nmax+1)
        for i = 0:nali.nmax
            w = get_transition_wvn(nali, i, nali.nmax)
            DOSp = get_DOS_matrix(E, L=26.5*1e4, g=20, νM=w + (10^(-6)), regime=2)
            #kplossp = BIRD.compute_kploss(Dfunc, x->BIRD.Dp(x, 26.5*1e4, w, 20), nali.nmax, rvals, e, C, 1000)
            #Jp = precision.(get_transport_matrix(E, V, DOSp, T, kploss=kplossp))
            Jp = precision.(get_transport_matrix(E, V, DOSp, T))
            kp = partialschur(Jp, nev=1, which=:SR)[1].eigenvalues[1]
            x[i+1] = kp/k0
        end
        lines!(ax, 0:nali.nmax, x, label="$T")
        scatter!(ax, 0:nali.nmax, x, label="$T")
    end
    
    axislegend(ax, L"T\;(K)", merge=true)
    fig
end