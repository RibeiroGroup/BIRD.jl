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

function get_transport_matrix(df, m::Morse{X}, T; kloss=1e8) where X
    E = transition_energy_matrix(m)
    V = transition_dipole_matrix(df, m)
    D = get_DOS_matrix(E)
    return get_transport_matrix(E, V, D, T; kloss=kloss)
end

function get_transport_matrix(E, V, D, T; kloss=1e8)
    
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
    end

    return J
end


