using Makie, WGLMakie, CairoMakie
using LaTeXStrings
using HDF5
using BIRD
using Interpolations
using Quadmath, DoubleFloats
using LinearAlgebra

### This is an auxiliary function to get E and V matrices (Transition energy and dipoles, repspectively)
### The default precision is Float128 from Quadmath. This is necessary given the small magnitude of the transitions
### at low temperature.
function get_EV(;precision=Float128)

    # Read data for dipole function
    path = joinpath(@__DIR__, "../QMcalcs/LiH/dip/dips.h5")
    r = vcat(0.0, h5read(path, "rvals"))
    dip = precision(0.529177) .* vcat(0.0, h5read(path, "dips")) # Convert from e⋅bohr to e⋅Å
    # Get interpolation
    itp = interpolate((r,), dip, Gridded(Linear()))
    # Dipole function. Return 0 for out of bounds
    lih_dipole(x) = x > 6.6 ? zero(precision) : itp(x)

    # Get Morse object for LiH
    lih = Morse(mA=precision(1.0)u"u", mB=precision(6.941)u"u", 
       ν=precision(1405.0)u"cm^-1", νχ=precision(23.1679)u"cm^-1", 
       re=precision(1.595)u"Å", De=precision(2.641)u"eV")

    # Get matrix with transition energies
    E = transition_energy_matrix(lih)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(lih_dipole, lih, rvals=0:0.01:30)

    return E, V
end

function hf_get_EV(;precision=Float128)

    # Dipole function
    hf_dipole(r) = 0.4541416928679838 * r * exp(-0.081616*r^4)

    # Get Morse object for LiH
    hf = Morse(mA=precision(1)u"u", mB=precision(19)u"u", 
        ν=precision(4138.0)u"cm^-1", νχ=precision(86.70)u"cm^-1", 
        re=precision(0.926)u"Å", De=precision(6.121646)u"eV")

    # Get matrix with transition energies
    E = transition_energy_matrix(hf)

    # Get matrix with transition dipole elements
    V = transition_dipole_matrix(hf_dipole, hf, rvals=0:0.01:30)

    return E, V
end

function temperature_analysis(E, V)
    fig = Figure(size=(500, 800), fontsize=22)
    ax1 = Axis(fig[1,1], ylabel=L"k_c/k_0", xticks=400:400:3000, 
    yticks=0.93:0.02:1.10, ylabelsize=30)

    ax2 = Axis(fig[2,1], xlabel="Temperature (K)", ylabel=L"k_p/k_0", xticks=400:400:3000, 
    yticks=0.0:30:300, ylabelsize=30)

    # Weak coupling
    #wc_temperature_analysis!(ax1, E, V, 9.6, 9.4)
    wc_temperature_analysis!(ax1, E, V, 3.31, 3.03)
    axislegend(ax1, merge=true, position=:rc)
    text!(ax1, 0.05, 0.98, text="(a)", align=(:left, :top), font=:bold, space=:relative)
    limits!(ax1, 1000, 3050, 0.92, 1.095)

    # Strong coupling
    #sc_temperature_analysis!(ax2, E, V, 9.6, 200, 1402, 1399)
    sc_temperature_analysis!(ax2, E, V, 3.31, 400, 3500, 3499)
    axislegend(ax2, merge=true, position=:rc)
    text!(ax2, 0.05, 0.98, text="(b)", align=(:left, :top), font=:bold, space=:relative)
    limits!(ax2, 1000, 3050, -5, 80)
    fig
end

function wc_temperature_analysis!(ax, E, V, L1, L2)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    DOSc1 = get_DOS_matrix(E, L=L1*1e4, regime=1)
    DOSc2 = get_DOS_matrix(E, L=L2*1e4, regime=1)

    Tvals = 1000:50:3000
    r1 = zeros(length(Tvals))
    r2 = zeros(length(Tvals))

    for i in eachindex(Tvals)
        T = Tvals[i]
        kc1 = eigvals(get_transport_matrix(E, V, DOSc1, T))[1]
        kc2 = eigvals(get_transport_matrix(E, V, DOSc2, T))[1]
        k0 = eigvals(get_transport_matrix(E, V, DOS0, T))[1]
        r1[i] = kc1 / k0
        r2[i] = kc2 / k0
        if T in [1000, 1500, 2000, 2500, 3000]
            println("T = $T  k0 = $(real(k0))   kc = $(real(kc1))")
        end
    end

    lines!(ax, Tvals, r1, label=L"$L_C = %$L1$ $\mu$m", color=Makie.wong_colors()[3])
    scatter!(ax, Tvals, r1, label=L"$L_C = %$L1$ $\mu$m", marker=:utriangle, color=Makie.wong_colors()[3])
    lines!(ax, Tvals, r2, label=L"$L_C = %$L2$ $\mu$m", color=Makie.wong_colors()[4])
    scatter!(ax, Tvals, r2, label=L"$L_C = %$L2$ $\mu$m", marker=:dtriangle, color=Makie.wong_colors()[4])
end

function sc_temperature_analysis!(ax, E, V, L, g, νM1, νM2)

    # Get free space density of states (regime = 0)
    DOS0 = get_DOS_matrix(E, regime=0)

    DOSp1 = get_DOS_matrix(E, L=L*1e4, g=g, νM=νM1, regime=2)
    DOSp2 = get_DOS_matrix(E, L=L*1e4, g=g, νM=νM2, regime=2)

    Tvals = 1000:50:3000
    r1 = zeros(length(Tvals))
    r2 = zeros(length(Tvals))

    for i in eachindex(Tvals)
        T = Tvals[i]
        kc1 = eigvals(get_transport_matrix(E, V, DOSp1, T))[1]
        kc2 = eigvals(get_transport_matrix(E, V, DOSp2, T))[1]
        k0 = eigvals(get_transport_matrix(E, V, DOS0, T))[1]
        if T in [1000, 1500, 2000, 2500, 3000]
            println("T = $T  k0 = $(real(k0))   kp = $(real(kc1))")
        end
        r1[i] = kc1 / k0
        r2[i] = kc2 / k0
    end

    lines!(ax, Tvals, r1, label=L"\omega_M = %$νM1")
    scatter!(ax, Tvals, r1, marker=:star5, label=L"\omega_M = %$νM1")
    lines!(ax, Tvals, r2, label=L"\omega_M = %$νM2")
    scatter!(ax, Tvals, r2, label=L"\omega_M = %$νM2")

    println("T = $(Tvals[1])  - kp ($νM1):  $(r1[1])   -  kp ($νM2):  $(r2[1])")
    println("T = $(Tvals[end])  - kp ($νM1):  $(r1[end])   -  kp ($νM2):  $(r2[end])")
end

function plank_at(ν, Tvals)

    fig = Figure()
    ax = Axis(fig[1,1])

    out = zeros(length(Tvals))
    ħ = BIRD.ħ
    k = BIRD.k
    ω = BIRD.ν2ω(ν)
    for i in eachindex(out)
        T = Tvals[i]
        out[i] = D0(ω) * ħ * ω * (1 / (exp(ħ*ω/(k*T)) - 1))
        out[i] = D0(ω) * ω * (1 / (exp(ħ*ω/(k*T)) - 1))
    end

    println(out)
    lines!(ax, Tvals, out)
    fig
end