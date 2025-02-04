using WGLMakie, Makie, CairoMakie
using LaTeXStrings
using HDF5
using Interpolations

function morse()

    fig = Figure(size=(600,800), fontsize=20)
    ax1 = Axis(fig[1,1], ylabel=L"Energy (eV)$$", title="(a) HF")
    ax2 = Axis(fig[2,1], xlabel=L"Bond Length (\AA)$$", ylabel=L"Energy (eV)$$", title="(b) LiH")

    hf = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV")

    lih = Morse(mA=1.0u"u", mB=6.941u"u", 
    ν=1405.0u"cm^-1", νχ=23.1679u"cm^-1", 
    re=1.595u"Å", De=2.641u"eV")

    # Morse energy
    E(m, r) = m.De * (1 - exp(-m.a*(r-m.re)))^2

    rvals = 0.0:0.01:20

    # HF
    hflevels = [0, 2, 5, 10, 15, 23]
    ν2eV = 	1.23981e-4
    for n in 0:hf.nmax
        En = ν2eV*energy(hf, n)
        rmax = hf.re - (1/hf.a) * log(1 - sqrt(En/hf.De))
        rmin = hf.re - (1/hf.a) * log(1 + sqrt(En/hf.De))
        lines!(ax1, [rmin, rmax < 4 ? rmax : 4.0], [En, En], color=(:lightblue4, 0.5))
    end

    for l in hflevels
        ψ = BIRD.compute_wfn(hf, l, rvals) 
        msk = abs2.(ψ) .> 1e-8 
        ψ .+= ν2eV*energy(hf, l)
        lines!(ax1, rvals[msk], ψ[msk], label=L"%$l")
    end
    lines!(ax1, rvals, [E(hf, r) for r in rvals], color=:midnightblue, linewidth=3)
    xlims!(ax1, 0, 7)
    ylims!(ax1, -0.2, 6.5)
    axislegend(ax1, L"$n$", position=:rb)


    ## LiH
    lihlevels = [0, 5, 10, 15, 20, 29]
    for n in 0:lih.nmax
        En = ν2eV*energy(lih, n)
        rmax = lih.re - (1/lih.a) * log(1 - sqrt(En/lih.De))
        rmin = lih.re - (1/lih.a) * log(1 + sqrt(En/lih.De))
        lines!(ax2, [rmin, rmax], [En, En], color=(:lightblue4, 0.5))
    end

    for l in lihlevels
        ψ = BIRD.compute_wfn(lih, l, rvals) 
        msk = abs2.(ψ) .> 1e-8 
        ψ .+= ν2eV*energy(lih, l)
        lines!(ax2, rvals[msk], ψ[msk], label=L"%$l")
    end
    lines!(ax2, rvals, [E(lih, r) for r in rvals], color=:midnightblue, linewidth=3)
    xlims!(ax2, 0, 12)
    ylims!(ax2, -0.2, 3.0)
    axislegend(ax2, L"$n$", position=:rb)

    fig
end

function dipoles()

    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4)

    path = joinpath(@__DIR__, "../QMcalcs/LiH/dip/dips.h5")
    r = vcat(0.0, h5read(path, "rvals")[4:end])
    dip = 0.529177 .* vcat(0.0, h5read(path, "dips")[4:end]) # Convert from e⋅bohr to e⋅Å
    # Get interpolation
    itp = interpolate((r,), dip, Gridded(Linear()))
    # Dipole function. Return 0 for out of bounds
    lih_dipole(x) = x > 6.6 ? 0.0 : itp(x)

    fig = Figure(fontsize=20)
    ax = Axis(fig[1,1], xlabel=L"Bond length (\AA)$$", ylabel=L"Dipole Moment ($e\cdot$\AA)",
    xticks=0:1:8)

    rvals = 0:0.05:10
    lines!(ax, 0.8:0.01:10, lih_dipole.(0.8:0.01:10), label="LiH (Interpolated)")
    scatter!(ax, r[2:end], dip[2:end], label="LiH (CCSD)")
    lines!(ax, rvals, hf_dipole.(rvals), label="HF")

    xlims!(ax, 0, 8)

    vlines!(ax, [0.926, 1.595], linestyle=:dash, color=:gray)

    text!(ax, 0.926, 0.7, text=L"r_e = 0.926\;\AA", align=(:center, :top), fontsize=18, rotation=-π/2)
    text!(ax, 1.595, 0.7, text=L"r_e = 1.595\;\AA", align=(:center, :bottom), fontsize=18, rotation=-π/2)

    axislegend(ax)

    return fig
end

function perturbation_analysis()

    fig = Figure(size=(650,400))
    gd = GridLayout(fig[1,1])
    axs = [Axis(gd[i,j]) for i = 1:2, j = 1:2]

    axs[1,1].xticks = 1:2:23
    axs[2,1].xticks = 1:2:23
    axs[1,1].yticks = 1:2:23
    axs[2,1].yticks = 1:2:23

    axs[1,2].xticks = 1:2:29
    axs[2,2].xticks = 1:2:29
    axs[1,2].yticks = 1:2:29
    axs[2,2].yticks = 1:2:29

    xlims!(axs[1,1], 10, 23.5)
    xlims!(axs[2,1], 10, 23.5)
    ylims!(axs[1,1], 10, 23.5)
    ylims!(axs[2,1], 10, 23.5)

    xlims!(axs[1,2], 16, 29.5)
    xlims!(axs[2,2], 16, 29.5)
    ylims!(axs[1,2], 16, 29.5)
    ylims!(axs[2,2], 16, 29.5)

    # Get molecular data
    hf = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV")

    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4)

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


    #### PLOT ANALYSIS FOR HF
    δ = 0.5 
    for i in 1:2

        m = i == 1 ? hf : lih
        dip = i == 1 ? hf_dipole : lih_dipole

        # Get matrix with transition energies
        E = transition_energy_matrix(m)

        # Get matrix with transition dipole elements
        V = transition_dipole_matrix(dip, m)

        # Get free space density of states (regime = 0)
        DOS = get_DOS_matrix(E, regime=0)

        for j in 1:2

            # Compute the sensitivity as S = k(after perturbation) / k(before perturbation) - 1
            δ = 1.5-(j-1)
            S = find_pertbation_matrix(E, V, DOS, T= i == 1 ? 4000 : 2000, δ=δ)

            #hm = heatmap!(axs[j,i], 0:m.nmax, 0:m.nmax, S, colormap=:bluesreds, colorrange=(-0.12, 0.12))
            hm = heatmap!(axs[j,i], 0:m.nmax, 0:m.nmax, S, colormap=:balance, colorrange=(-0.12, 0.12))
            if i == 2 && j == 2
                Colorbar(gd[:, end+1], hm, label="Sensitivity")
            end
            #translate!(hm, 0, 0, -100)

            mol = i == 1 ? "HF" : "LiH"
            lt = '`' + i + 2*(j÷2)
            text!(axs[j,i], 0.02, 0.9, text=L"(%$lt) %$mol, $\delta = $ %$δ", align=(:left, :top), space=:relative)

        end
    end

    Label(gd[1:2,0], L"i", fontsize=25)
    Label(gd[3,1:2], L"j", fontsize=25)

    rowgap!(gd, 2, 1.5)

    fig
end

function perturbation_table(mol; overtones=true)

    # Get molecular data
    hf = Morse(mA=1u"u", mB=19u"u", 
    ν=4138.0u"cm^-1", νχ=86.70u"cm^-1", 
    re=0.926u"Å", De=6.121646u"eV")

    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4)

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

    m = lowercase(mol) == "hf" ? hf : lih
    dipf = lowercase(mol) == "hf" ? hf_dipole : lih_dipole
    # Get matrix with transition energies
    E = transition_energy_matrix(m)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(dipf, m, overtones=overtones)

    # Get free space density of states (regime = 0)
    DOS = get_DOS_matrix(E, regime=0)


    out = """\\begin{table}[h]
    \\centering
    \\caption{Most important transitions according to their sensitivity analysis score for $(mol).}
    \\label{tab:$(mol)_sensitivity}
    \\begin{tabular}{|c|c|c|c|}
    \\hline
    \\textbf{Transition} & \\textbf{Energy (cm\$^{-1}\$)} & \\textbf{Sensitivity} (\$\\delta = 1.5\$) & \\textbf{Sensitivity} (\$\\delta = 0.5\$) \\\\
    \\hline
    \\hline"""

    S1 = find_pertbation_matrix(E, V, DOS, T= lowercase(mol) == "hf" ? 4000 : 2000, δ=1.5)
    S2 = find_pertbation_matrix(E, V, DOS, T= lowercase(mol) == "hf" ? 4000 : 2000, δ=0.5)
    l1, Sl1 = _linearize_offdiagonal(S1)
    l2, Sl2 = _linearize_offdiagonal(S2)
    pe1 = reverse(sortperm(abs.(Sl1)))
    pe2 = reverse(sortperm(abs.(Sl2)))
    for k = 1:10
        #@assert l1[pe1][k] == l2[pe2][k]
        (i,j) = l1[pe1][k]
        s1 = round(Sl1[pe1][k], digits=3)
        s2 = round(Sl2[pe1][k], digits=3)
        wvn = round(BIRD.get_transition_wvn(m, i, j), digits=1)

        new = "\$$i \\rightarrow $j\$ & $wvn & $s1 & $s2 \\\\"
        out *= "\n"
        out *= new
        out *= "\n"
        out *= "\\hline"
    end

    out *= """

    \\end{tabular}
    \\end{table}"""

    print(out)
end


function find_pertbation_matrix(E, V, DOS; δ=1.2, T=300)

    # Get free space rate
    J0 = get_transport_matrix(E, V, DOS, T)
    k0 = eigmin(J0)

    N = size(E,1)

    S = similar(E)
    S .= 0.0

    for i in 1:N
        S[i,i] = 0.0
        for j in (i+1):N
            DOS[i,j] *= δ
            DOS[j,i] *= δ
            J = get_transport_matrix(E, V, DOS, T)
            S[i,j] = eigmin(J) / k0 - 1
            S[j,i] = S[i,j]
            DOS[i,j] /= δ
            DOS[j,i] /= δ
        end
    end

    return S
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