export D0, Dc, Dp, get_DOS_matrix

"""
    get_DOS_matrix(E; L=20000, νM=3000, g=200, regime=0)

Given a symmetric matrix E with transition frenquencies (angualr, in Hz), computes a corresponding
matrix with density of states at each particular transition energy, in units of Hz⁻¹Å⁻³.

The keyword `regime` commands which type of density of states is computed:

| regime | Type of DOS                     |
|--------|---------------------------------|
| 0      | Free-space                      |
| 1      | Planar Cavity (weak coupling)   |
| 2      | Planar Cavity (Strong coupling) |

Keyword arguments L, νM, g are necessary depending on the `regime`:

| Keyword | Description                                   | Std. Value | Units | Regime |
|---------|-----------------------------------------------|------------|-------|--------|
| L       | Distance between mirrors in the planar cavity | 10000      | Å     | 1,2    |
| νM      | Strongly coupled transition energy            | 3000       | cm⁻¹  | 2      |
| g       | Light-matter coupling strength                | 200        | cm⁻¹  | 2      |

These keywords are unused if the corresponding `regime` does not require it. 
"""
function get_DOS_matrix(E::Matrix{T}; L=20000, νM=3000, g=200, regime=0) where T

    N = size(E, 1)
    DOS = zeros(T, N, N)

    f(ω, L, νM, g) = if regime == 0
        return D0(ω)
    elseif regime == 1
        return Dc(ω, L)
    elseif regime == 2
        return Dp(ω, L, ν2ω(νM), ν2ω(g))
    end

    for i in axes(DOS, 1)
        for j in (i+1):N
            DOS[i,j] = f(E[i,j], L, νM, g)
            DOS[j,i] = DOS[i,j]
        end
    end
        
    return DOS
end

"""
    D0(ω)

Returns the EM density of states in free spaces at the angular frequency ω. 
Input ω must be in rad.s⁻¹ and the output is given in Hz⁻¹Å⁻³.
"""
function D0(ω)
    ω^2 / (π^2 * c^3)
end

"""
    D0(ω)

Returns the EM density of states in free spaces at the angular frequency ω. 
Input ω must be in rad.s⁻¹ and L in Å and the output is given in Hz⁻¹Å⁻³.
"""
function Dc(ω, L)
    ω / (π * c^2 * L) * ( floor(ω*L/(π*c)) + 0.5 )
end

"""
    Dp(ω, L, ωM, g)

Computes the photonic weighted polariton density of states at the angular frequency ω.
L is the cavity length in Å, ωM is the matter angular frequency (strongly coupled to the EM field),
and g is the coupling strength also in units of angular frequency (rad/s).
The output is given in Hz⁻¹Å⁻³.
"""
function Dp(ω, L, ωM, g; verbose=false)

    if ω < ωM && (ωM/ω - 1) < 1e-10
        @warn "Density of states at ω = ωM goes to ±∞ in the polaritonic case."
        return Inf
    end

    # Print function for debbuging
    printif(t) = verbose ? println(t) : nothing

    printif(repeat('-', 30))
    printif("• Input")
    printif("Probing Frequency:")
    printif("ω: $(round(ω*1e-12, digits=2)) THz")
    printif("ν: $(round(ω2ν(ω), digits=2)) cm⁻¹")
    printif("\nMatter Frequency:")
    printif("ωM: $(round(ωM*1e-12, digits=2)) THz")
    printif("νM: $(round(ω2ν(ωM), digits=2)) cm⁻¹")
    printif("\nCoupling Strength:")
    printif("g: $(round(g*1e-12, digits=2)) THz")
    printif("νM: $(round(ω2ν(g), digits=2)) cm⁻¹")
    printif("\nCavity:")
    printif("Length: $(L*1e-8) cm - $(L*1e-4) μm - $(L) Å")
    printif("Cutoff Energy: $(round(c*π/L*1e-12, digits=2)) THz - $(round(1/(2*L*1e-8), digits=2)) cm⁻¹")
    printif(repeat('-', 30))

    # Convert L from Å to cm
    L = L * 1e-8

    # Compute cavity q0 in cm
    m = 0
    q0 = m*π/L

    # For the given frequency and m value, find the matching cavity energy
    ωc2 = ω^2 * (ωM^2 + g^2 - ω^2) / (ωM^2 - ω^2)

    # ωc2 will be negative if ω > ωM AND ω < √ωM + g. In that case returns zero. (i.e. there is no valid q)
    if ωc2 < 0.0
        SGl = "$(round(ω2ν(ωM)))"
        SG = "$(round(ω2ν(sqrt(ωM^2 + g^2))))"
        printif("Probed frequency ($(round(ω2ν(ω),digits=2)) cm⁻¹ is inside the stop gap region ($SGl - $SG ) cm⁻¹")
        return 0.0
    end

    ωc = sqrt(ωc2)

    printif("\nCavity Energy for the input frequency:")
    printif("ωc: $(round(ωc*1e-12, digits=2)) THz")
    printif("νc: $(round(ω2ν(ωc), digits=2)) cm⁻¹")

    out = 0.0

    printif("------ Branch contributions -----")
    # Condition: cavity cutoff energy can't be larger than the desired ω
    while c_cm*q0 < ωc

        # Reverse the equation for ωc to obtain q in cm⁻¹
        q = sqrt( ωc2/c_cm^2 - q0^2 )

        # Compute photonic content of the mode
        Ω = (ω^2 - ωM^2 - g^2)^2
        Pc = Ω / (Ω + ωc2*g^2)

        # Compute group velocity: This will depend on whether the mode is LP or UP
        # We will considet that ωLP < ωM and ωUP > ωM
        # vg has units of cm/s
        if ω > ωM
           vg = q*c_cm^2/(2*ω) * (1 +  (ωc^2 - ωM^2 + g^2)/(sqrt( (ωc^2 + ωM^2 +g^2)^2 - 4ωc^2 * ωM^2)))
        elseif ω < ωM
           vg = q*c_cm^2/(2*ω) * (1 -  (ωc^2 - ωM^2 + g^2)/(sqrt( (ωc^2 + ωM^2 +g^2)^2 - 4ωc^2 * ωM^2)))
        else
            @warn "Density of states at ω = ωM goes to ±∞ in the polaritonic case."
            return Inf
        end

        # Put results together: q * Pc / π*L*vg
        branch = (m == 0 ? 0.5 : 1) * q*Pc / (π*L*vg)
        out += branch

        printif("m = $m | q₀ = $(round(ω2ν(c_cm*q0), digits=2)) cm⁻¹")

        # Update m and q0
        m += 1
        q0 = m*π/L

    end
    printif("---------------------------------")

    printif("Number of bands contributing to the DOS at ω = $(round(ω*1e-12, digits=2)) THz (ν = $(round(ω2ν(ω),digits=2)) cm⁻¹) = $(m)")
    # Convert output from s/cm^3 to s/Å^3
    return out * 1e-24
end