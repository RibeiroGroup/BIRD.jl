using BIRD
using DataInterpolations
using LinearAlgebra
using Quadmath, DoubleFloats, ArnoldiMethod

function perturbation_analysis()

    fig = Figure(size=(650,400))
    gd = GridLayout(fig[1,1])
    ax = Axis(gd[1,1], xlabel=L"j", ylabel=L"i")

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

    # Compute the sensitivity as S = k(after perturbation) / k(before perturbation) - 1
    S = find_pertbation_matrix(E, V, DOS, T=2000, δ=1.5)

    hm = heatmap!(ax, 0:nali.nmax, 0:nali.nmax, S, colormap=:balance, colorrange=(-0.12, 0.12))
    Colorbar(gd[:, end+1], hm, label="Sensitivity")

    fig
end

function overtone_sensitivity(;T=2000, δ=1.5, precision=Float64)

    fig = Figure(size=(650,400))
    ax = Axis(fig[1,1], xlabel=L"i", ylabel="Sensitivity")

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

    # Get free space rate
    J0 = precision.(get_transport_matrix(E, V, DOS, T))
    k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
    # k0 = eigmin(J0)

    j = nali.nmax
    S = zeros(j)

    for i in 1:j
        DOS[i,j+1] *= δ
        DOS[j+1,i] *= δ
        J = precision.(get_transport_matrix(E, V, DOS, T))
        kp = partialschur(J, nev=1, which=:SR)[1].eigenvalues[1]
        # kp = eigmin(J)
        S[i] = kp / k0 - 1
        DOS[i,j+1] /= δ
        DOS[j+1,i] /= δ
    end

    println(j)
    println(size(DOS))

    scatter!(ax, 0:(j-1), S)
    fig
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