using BIRD
using LinearAlgebra
using DataInterpolations
using Quadmath, DoubleFloats, ArnoldiMethod
using HDF5

function lowest_eigenval(M::Matrix{Float64})
    return eigmin(M)
end

function lowest_eigenval(M::Matrix{Float128})
    return partialschur(M, nev=1, which=:SR)[1].eigenvalues[1]
end

function generate_strong_coupling(outname; νMvals=200:0.5:1000, L=5, T=400, precision=Float128, g=125)

    #### Get NaLi data
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "../../Fedorov2014.txt"))]
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
        DOSp = get_DOS_matrix(E, L=L*μm2Å, regime=1) # NOTE THAT I CHANGED IT TO THE WEAK COUPLING CONDITION TO ISOLATE THE EFFECT OF KPLOSS
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
