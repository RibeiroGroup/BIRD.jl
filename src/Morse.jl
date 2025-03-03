using LinearAlgebra

export Morse, energy, get_fundamentals, get_transition_angfrequency, get_transition_frequency, get_transition_wvn
export get_dipole, compute_wfn, transition_dipole, transition_dipole_matrix, transition_energy_matrix, potential

struct Morse{T}
    redmass::T
    ν::T
    νχ::T
    re::T
    De::T
    a::T
    nmax::Int
    λ::T
end

"""
    mA and mB are the masses of the atoms A and B in a.m.u
"""
function Morse(;mA, mB, ν, νχ, re, De)

    c = CODATA2018.SpeedOfLightInVacuum
    h = CODATA2018.PlanckConstant
    e = CODATA2018.ElementaryCharge
    ħ = h / 2π

    # Reduced mass (in a.m.u)
    redmass = mA * mB / (mA + mB)

    # Force constant using the relation ν = (1/2π) * √(ke/μ) where ν is the frequency
    ke = redmass * (2π * c * ν)^2 

    # Morse parameters a from force constant and dissociation energy
    a = sqrt(ke / (2*De)) 

    # Compute number of bound states
    λ1 = 2*redmass*De 
    λ2 = ħ*a
    λ = √λ1 / λ2 |> upreferred

    nmax = Int(floor(λ - 0.5))

    return Morse(ustrip(u"u", redmass), 
                ustrip(u"cm^-1", ν),
                ustrip(u"cm^-1", νχ),
                ustrip(u"Å", re), 
                ustrip(u"eV", De), 
                ustrip(u"Å^-1", a), nmax, λ)
end

"""
# energy(M::Morse{T}, n)
    Returns the n-th level energy (in cm⁻¹) of the morse oscillator
"""
function energy(M::Morse{T}, n) where T
    return M.ν * (n + 0.5) - M.νχ * (n + 0.5)^2
end

"""
# potential(M::Morse{T}, r)
    Computes the Morse potential energy (in cm⁻¹) at the displacement r (in Å)
"""
function potential(M::Morse{T}, r) where T
    return eV2cm * M.De * (1 - exp(-M.a*(r-M.re)))^2
end

"""
# get_transition_wvn(M::Morse{T}, i, j)
    Computes the transition energy (in wavenumbers) i -> j
"""
function get_transition_wvn(M::Morse{T}, i, j) where T
    return abs(energy(M, j) - energy(M, i))
end

"""
# get_transition_frequency(M::Morse{T}, i, j)
    Computes the transition energy (in frequency) i -> j
"""
function get_transition_frequency(M::Morse{T}, i, j) where T
    c = ustrip(u"cm/s", CODATA2018.SpeedOfLightInVacuum)
    return c * get_transition_wvn(M, i, j)
end

"""
# get_transition_angfrequency(M::Morse{T}, i, j)
    Computes the transition energy (in angular frenquency) i -> j
"""
function get_transition_angfrequency(M::Morse{T}, i, j) where T
    return 2π*get_transition_frequency(M, i, j)
end

"""
# compute_wfn(M::Morse{T}, n, rvals)

    Compute the Morse wavefunction for the quantum number `n` over a range of values `rvals`.

The evaluation uses the formula found [here.](https://en.wikipedia.org/wiki/Morse_potential)
"""
function compute_wfn(M::Morse{T}, n, rvals) where T

    out = zeros(T, length(rvals))

    for i in eachindex(rvals)
        r = rvals[i]
        z = 2*M.λ * exp(-M.a*(r - M.re))
        α = 2*M.λ - 2*n - 1

        out[i] = z^(M.λ - n - 0.5) * exp(-0.5*z) * laguerre(α, n, z)
    end

    normalize!(out)
    if abs(out[end]) > 1e-6
        @warn "Ψ = ($(out[end])) > 1e-6 at the end of the chosen interval! This may raise integration errors"
    end
    if abs(out[1]) > 1e-6
        @warn "Ψ = ($(out[1])) > 1e-6 at the beginning of the chosen interval! This may raise integration errors"
    end

    return out
end

"""
# laguerre(α, n, x)

    Evaluate the generalized Laguerre polynomial L[α,n] at the point x. 
The evaluation is made using the recursive formula found [here.](https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials)
"""
function laguerre(a, n, x)
    if n == 0
        return 1.0
    elseif n == 1
        return 1.0 + a - x
    else
        L0 = 1.0
        L1 = 1.0 + a - x
        for k in 2:n
            L = ((2k - 1 + a - x) * L1 - (k - 1 + a) * L0) / k
            L0, L1 = L1, L
        end
        return L1
    end
end

"""
    transition_dipole(df, M, i, j)
    transition_dipole(df, M, i, j, rvals)

Computes the transition dipole moment ⟨i|μ|j⟩ for the Morse system `M` using the dipole
function `df` with `rvals` as the integration interval for the morse function. 
If `rvals` is omitted, the function `find_rvals` will be used to determine the minimal range.
"""
function transition_dipole(df, M::Morse{T}, i, j, rvals) where T

    ψ1 = compute_wfn(M, i, rvals)
    ψ2 = compute_wfn(M, j, rvals)

    for i in eachindex(rvals)
        ψ2[i] *= df(rvals[i])
    end

    return dot(ψ1, ψ2)
end

function transition_dipole(df, M::Morse{T}, i, j) where T
    # Find appropriate range for integration
    rvals = find_rvals(M, max(i,j))
    transition_dipole(df, M, i, j, rvals)
end

"""
    find_rvals(M, n)

Find the minimal integration interval for the morse function with quantum number n.
The criteria used is that the value of the wave function must be < 1e-6 at the integration boundaries.
"""
function find_rvals(M::Morse{T}, n) where T

    # Get rmin (0.0 or 0.5)
    rmax = 3.0

    z = 2*M.λ * exp(-M.a*(rmax - M.re))
    α = 2*M.λ - 2*n - 1

    wfn = z^(M.λ - n - 0.5) * exp(-0.5*z) * laguerre(α, n, z)

    while abs(wfn) > 1e-6

        if rmax > 30.0
            error("Could not find a upper bound for integration for the state n = $n (rmax > 30)")
        end

        rmax += 1
        z = 2*M.λ * exp(-M.a*(rmax - M.re))
        α = 2*M.λ - 2*n - 1
        wfn = z^(M.λ - n - 0.5) * exp(-0.5*z) * laguerre(α, n, z)
    end

    return 0.0:0.01:rmax
end

function transition_dipole_matrix(df, M::Morse{T}; overtones=true, rvals=nothing) where T

    V = zeros(T, M.nmax+1, M.nmax+1)

    findr = isnothing(rvals)

    if overtones
        for i = 1:(M.nmax+1)
            for j = i:(M.nmax+1)
                V[i,j] = findr ? transition_dipole(df, M, i-1, j-1) : transition_dipole(df, M, i-1, j-1, rvals)
                V[j,i] = V[i,j]
            end
        end

        return V
    else
        for i = 1:(M.nmax)
            j = i+1
            V[i,j] = findr ? transition_dipole(df, M, i-1, j-1) : transition_dipole(df, M, i-1, j-1, rvals)
            V[j,i] = V[i,j]
        end

        return V
    end
end

function transition_energy_matrix(M::Morse{T}) where T
    E = zeros(T, M.nmax+1, M.nmax+1)

    for i = 1:(M.nmax+1)
        for j = i:(M.nmax+1)
            E[i,j] = get_transition_angfrequency(M, i-1, j-1)
            E[j,i] = E[i,j]
        end
    end

    return E
end

function get_continuum_states(m::Morse; N=2000, rmin=0.0, rmax=30.0)

    # Potential energy function for the numerov procedure, in cm⁻¹
    Vf(x) = potential(m, x)

    # Prefactor so units work out - output of numerov will be cm⁻¹
    pf = get_wvn_prefactor(u"Å", m.redmass * PhysicalConstants.CODATA2018.AtomicMassConstant)

    # Get energies and wavefunctions from numerov
    println("Running Numerov...")
    e, C = numerov(Vf, N, rmin, rmax, pf)
    println("Done!")

    # Compute dipole over the grid
    rvals = range(start=rmin, stop=rmax, length=N)

    return rvals, e, C
end