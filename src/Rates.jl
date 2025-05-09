export get_transport_matrix

function spontaneous_emission_rate(D, ω, μ)
    ω * μ^2 * π / (3*ε0*ħ) * D
end

function absorption_rate(D, T, ω, μ)
    # kabs = kstim
    ω * μ^2 * π / (3*ε0*ħ) * D * (1 / (exp(ħ*ω/(k*T)) - 1))
end

function stimulated_emission_rate(D, T, ω, μ)
    ω * μ^2 * π / (3*ε0*ħ) * D * (1 / (exp(ħ*ω/(k*T)) - 1))
end

function get_transport_matrix(df, m::Morse{X}, T; kloss=1e8, kploss=false) where X
    E = transition_energy_matrix(m)
    V = transition_dipole_matrix(df, m)
    D = get_DOS_matrix(E)
    return get_transport_matrix(E, V, D, T; kloss=kloss, kploss=kploss)
end

function get_transport_matrix(E, V, D, T; kloss=1e8, kploss=false)
    
    N = size(E,1)
    # Intialize the transport matrix
    J = similar(E) 
    J .= 0.0

    for i = axes(E,1)
        # Get absorption terms j -> i
        for j = 1:(i-1)
            # Note that we add 1 to the indexes because while the states are labelled from 0, the arrays start at 1
            J[i, j] -= absorption_rate(D[i,j], T, E[i,j], V[i,j]) 
        end

        # Get emission terms j -> i
        for j = (i+1):N
            J[i, j] -= spontaneous_emission_rate(D[i,j], E[i,j], V[i,j])
            J[i, j] -= stimulated_emission_rate(D[i,j], T, E[i,j], V[i,j])
        end

        # Get diagonal terms i -> k
        for k in 1:(i-1)
            J[i, i] += spontaneous_emission_rate(D[i,k], E[i,k], V[i,k])
            J[i, i] += stimulated_emission_rate(D[i,k], T, E[i,k], V[i,k])
        end

        for k in (i+1):N
            J[i, i] += absorption_rate(D[i,k], T, E[i,k], V[i,k]) 
        end

        if i == N
            J[i, i] += kloss
        end

        # Add photodissociation to continuum
        if kploss != false
            J[i,i] += kploss[i]
        end
    end

    return J
end

function compute_kploss(dipole_function, DOS_function, m::Morse, T; N=2000, rmin=0.0, rmax=30.0)
    rvals, e, C = get_continuum_states(m, N=N, rmin=rmin, rmax=rmax)
    compute_kploss(dipole_function, DOS_function, m.nmax, rvals, e, C, T)
end

function compute_kploss(dipole_function, DOS_function, nmax, rvals, e, C, T)
    μAi = continuum_transition_dipoles(dipole_function, nmax, rvals, e, C)
    compute_kploss(DOS_function, e, μAi, T)
end

function continuum_transition_dipoles(dipole_function, nmax, rvals, e, C)

    # Compute dipole moments over the grid values
    dips = dipole_function.(rvals)

    # Allocate intermediate array
    ψxd = zeros(length(rvals))

    # Number of bound states
    Nbound = nmax + 1

    # Number of continuum states
    Ncontinuum = length(e) - Nbound

    # Allocate output Matrix
    μAi = zeros(Ncontinuum, Nbound)

    # Loop over bound states - i
    for i in axes(μAi, 2)
        ψxd .= C[:,i] .* dips

        # Loop over continuum states - A
        for A in axes(μAi, 1)
            # Compute index of the continuum state
            j = A + Nbound
            # Compute transition dipole moment 
            μAi[A,i] = dot(ψxd, C[:,j])
        end
    end

    return μAi
end

function compute_kploss(DOS_function, e, μAi, T)

    # Allocate output array
    phot_losses = zeros(size(μAi, 2))

    # Loop over bound states
    for i in axes(μAi, 2)

        # Loop over continuum states
        for A in axes(μAi, 1)

            # Compute index of the continuum state
            j = A + size(μAi, 2)
            # Compute transition angular frequency
            ω = ν2ω(e[j] - e[i])
            # Density of state at that frequency
            D = DOS_function(ω)

            # Compute and accumulate absorption rates (will sum over all unbounded states)
            phot_losses[i] += absorption_rate(D, T, ω, μAi[A,i])
        end
    end

    return phot_losses
end