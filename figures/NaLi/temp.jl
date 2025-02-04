using LinearAlgebra

include("get_pes.jl")

function nali_kps_vs_T()

    fig = Figure()
    ax = Axis(fig[1,1])

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

    μm2Å = 1e4

    DOS0 = get_DOS_matrix(E, regime=0)
    DOSc = get_DOS_matrix(E, L=9.6*μm2Å,  regime=1)
    DOSp = get_DOS_matrix(E, L=9.6*μm2Å, ωM=1.0,regime=2)


end

function lowT_freespace_limit()
    fig = Figure()
    ax1 = Axis(fig[1,1])
    #ax2 = Axis(fig[1,2])

    three_level = Morse(mA=22.98u"u", mB=6.941u"u", 
           ν=750.0u"cm^-1", νχ=39.245u"cm^-1", 
           re=2.889u"Å", De=0.16u"eV")

    # Dipole function for HF
    dipolef(r) = 0.4541416928679838 * r * exp(-0.081616*r^4)
    E = transition_energy_matrix(three_level)
    V = transition_dipole_matrix(dipolef, three_level, rvals=0.0:0.01:50, overtones=true)

    DOS0 = get_DOS_matrix(E, regime=0)

    Tvals = 70:10:400

    A02s = zeros(length(Tvals))
    k0s = zeros(length(Tvals))
    kps = zeros(length(Tvals))
    approx = zeros(length(Tvals))

    for i in eachindex(Tvals)
        T = Tvals[i]
        J = get_transport_matrix(E, V, DOS0, T)
        k0s[i] = eigmin(J)
        A02s[i] = -J[3,1]

        α = 1000
        Ja = alpha_Jp(α, E, V, DOS0, T)
        kps[i] = eigmin(Ja)

        A01 = -Ja[2,1]
        A02 = -Ja[3,1]
        A12 = -Ja[3,2] / α
        E10 = -Ja[1,2]

        approx[i] = (α*(A01*A12 + A02*A12) + E10*A02) / (E10 + α*A12)
    end

    scatter!(ax1, Tvals, k0s, label="Numerical Diagonalization α = 1.0")
    lines!(ax1, Tvals, A02s, label="Approximation as A02")

    scatter!(ax1, Tvals, kps, label="Numerical Diagonalization α = 10³")
    lines!(ax1, Tvals, approx, label="Approximation")
    #axislegend(ax2)
    axislegend(ax1, position=:lt)

    fig
end


function analytic_eigen(E, V, DOS0, T, α; lowT=false, lowa=false, higha=false)
    J = alpha_Jp(α, E, V, DOS0, T)
    analytic_eigen(J, α, lowT=lowT, lowa=lowa, higha=higha)
end

function analytic_eigen(J, α; lowT=false, lowa=false, higha=false)
    A01 = -J[2,1]
    A02 = -J[3,1]
    A12 = -J[3,2] / α
    E10 = -J[1,2]
    E20 = -J[1,3]
    E21 = -J[2,3] / α

    a = 1
    b = - (A01 + A02 + E10 + α*A12)
    c = (A01 + A02)*(E10 + α*A12) - E10*A01
    Δ = b^2 - 4*a*c
    println("b = $b")
    println("c = $c")
    println("Δ = $Δ")

    if lowT && lowa
        srs = A02 #+ A02^2 / E10 
        println("Low T low a Series approx: $srs")
    end

    if lowT && higha
        N = α * (A01*A12 + A02*A12) + A02*E10
        D = E10 + α*A12
        println("\nLow T high a Series approx: $(N/D)")
    end

    λ1 = (-b + sqrt(Δ)) / (2*a)
    λ2 = (-b - sqrt(Δ)) / (2*a)

    return min(λ1, λ2)
end

function alpha_Jp(alpha, E, V, DOS0, T)
    J = get_transport_matrix(E, V, DOS0, T)
        
    J[2,2] = -J[1,2] -alpha*J[3,2]
    J[3,2] = alpha*J[3,2]
    J[2,3] = alpha*J[2,3]
    return J
end

function lowTrates(J)
    A01 = -J[2,1]
    A02 = -J[3,1]
    A12 = -J[3,2] 
    E10 = -J[1,2]

    N = (A01*A12 + A02*A12) + A02*E10
    D = E10 + A12

    return N/D
end

function ratio(a1, a2, E, V, DOS0, T)
    J = get_transport_matrix(E, V, DOS0, T)
    A01 = -J[2,1]
    A02 = -J[3,1]
    A12 = -J[3,2] 
    E10 = -J[1,2]

    N = (a1*(A01*A12 + A02*A12) + A02*E10) * (E10 + a2*A12)
    D = (a2*(A01*A12 + A02*A12) + A02*E10) * (E10 + a1*A12)

    return N/D
end

function approx_ratio(a1, a2, E, V, DOS0, T)
    J = get_transport_matrix(E, V, DOS0, T)
    A01 = -J[2,1]
    A02 = -J[3,1]
    A12 = -J[3,2] 
    E10 = -J[1,2]

    N = (a1*(A01*A12 + A02*A12) + A02*E10)
    D = (A02) * (E10 + a1*A12)

    return N/D
end

function get_exact_k(a, Tvals)
    ks = [eigmin(alpha_Jp(a, E, V, DOS0, T)) for T in Tvals]
    return ks
end

function plot_ratios(a1, a2, Tvals)
    k1 = get_exact_k(a1, Tvals)
    k2 = get_exact_k(a2, Tvals)
    Plots.scatter(Tvals, k1 ./ k2, label="Exact")
end

function plot_ratios!(a1, a2, Tvals)
    k1 = get_exact_k(a1, Tvals)
    k2 = get_exact_k(a2, Tvals)
    Plots.scatter!(Tvals, k1 ./ k2, label="Exact")
end

function plot_approx_ratios!(a1, a2, Tvals)
    r = [ratio(a1, a2, E, V, DOS0, T) for T in Tvals]
    Plots.scatter!(Tvals, r, label="Approx")
end

function plot_approx_ratios(a1, a2, Tvals)
    r = [ratio(a1, a2, E, V, DOS0, T) for T in Tvals]
    r2 = [approx_ratio(a1, a2, E, V, DOS0, T) for T in Tvals]
    Plots.scatter(Tvals, r, label="Full")
    Plots.scatter!(Tvals, r2, label="Approx")
end

function plot_elements(E, V, DOS0, Tvals)

    A01s = zeros(length(Tvals))
    A12s = zeros(length(Tvals))
    A02s = zeros(length(Tvals))
    E10s = zeros(length(Tvals))

    for i in eachindex(Tvals)
        J = get_transport_matrix(E, V, DOS0, Tvals[i])
        A01s[i] = -J[2,1]
        A02s[i] = -J[3,1]
        A12s[i] = -J[3,2] 
        E10s[i] = -J[1,2]
    end

    Plots.scatter(Tvals, A01s, label="A01")
    Plots.scatter!(Tvals, A12s, label="A12")
    Plots.scatter!(Tvals, A02s, label="A02")
    #Plots.scatter!(Tvals, E10s, label="E10")
end


function BEfactor(e, T)
    1 / (exp(BIRD.ħ*e/(BIRD.k*T)) - 1)
end

function BEapprox1(e, T; order=2)

    a = BIRD.ħ*e/BIRD.k
    out = a/T

    for n in 2:order
        out += a^n / (factorial(n) * T^n)
    end

    return 1/out
end

function BEapprox2(e, T; order=2)
    exp(-BIRD.ħ*e/(BIRD.k*T))
end

function ratio_approx(α, ϵ01, ϵ12, Tvals)
    rout = zeros(length(Tvals))

    for i in eachindex(rout)
        β = BIRD.ħ/(BIRD.k*Tvals[i])
        rout[i] = 1 + α / (1 + exp(-β*ϵ01) + α*exp(-β*ϵ12))
    end

    return rout
end

function plot_anal_ratios(αvals, ϵ01, ϵ12, Tvals)

    Plots.plot()
    for α in αvals
        r = ratio_approx(α, ϵ01, ϵ12, Tvals)

        Plots.plot!(Tvals, r, label="α = $α", linewidth=5)
    end
    Plots.plot!()
end
