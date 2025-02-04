using LinearAlgebra
using LaTeXStrings

# Gets a 3x3 transport matrix where one transition is scaled by α
function alpha_Jp(alpha, E, V, DOS0, T)

    @assert size(E) == size(V) == size(DOS0) == (3,3)

    J = get_transport_matrix(E, V, DOS0, T)
        
    J[2,2] = -J[1,2] -alpha*J[3,2]
    J[3,2] = alpha*J[3,2]
    J[2,3] = alpha*J[2,3]
    return J
end

function lowT_rates()
    fig = Figure(fontsize=20)
    ax1 = Axis(fig[1,1], xlabel=L"Temperature (K)$$", ylabel=L"Rate (s$^{-1}$)")

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

    scatter!(ax1, Tvals, k0s, label=L"Numerical Diagonalization $α = 1.0$")
    lines!(ax1, Tvals, A02s, label=L"Approximation $\lambda_1$")

    scatter!(ax1, Tvals, kps, label=L"Numerical Diagonalization $α = 10^3$")
    lines!(ax1, Tvals, approx, label=L"Approximation $\lambda_\alpha$")
    axislegend(ax1, position=:lt)

    fig
end

function three_level_temp_analysis()

    fig = Figure(fontsize=20, size=(800,400))
    ax1 = Axis(fig[1,1], ylabel=L"Ratio $r = \lambda_\alpha / \lambda_1$", title="(a) Approximation")
    ax2 = Axis(fig[1,2], title="(b) Numerical Diagonalization")
    Label(fig[2,1:2], L"Temperature (K)$$")

    three_level = Morse(mA=22.98u"u", mB=6.941u"u", 
           ν=750.0u"cm^-1", νχ=39.245u"cm^-1", 
           re=2.889u"Å", De=0.16u"eV")

    # Dipole function for HF
    dipolef(r) = 0.4541416928679838 * r * exp(-0.081616*r^4)
    E = transition_energy_matrix(three_level)
    V = transition_dipole_matrix(dipolef, three_level, rvals=0.0:0.01:50, overtones=true)
    DOS0 = get_DOS_matrix(E, regime=0)

    ϵ = BIRD.ħ*(E[1,2] + E[2,3])/2 # units: J

    Tvals = 70:10:400


    for α in [1, 10, 100, 1000]
        rs = zeros(length(Tvals))
        r_ex = zeros(length(Tvals))
        for i in eachindex(Tvals)
            T = Tvals[i]
            β = 1/(BIRD.k*T) # units: J
            rs[i] = 1 + α / (2 * (1 + exp(-β*ϵ)*(1+2α) ))


            # Exact calc
            J = get_transport_matrix(E, V, DOS0, T)
            Jα = alpha_Jp(α, E, V, DOS0, T)
            r_ex[i] = eigmin(Jα) / eigmin(J)
        end
        lines!(ax1, Tvals, rs, linewidth=2, label=L"%$α")
        lines!(ax2, Tvals, r_ex, linewidth=2, label=L"%$α")
    end

    Legend(fig[1,3], ax1, L"\alpha")
    fig

end