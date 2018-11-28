#!/usr/bin/env julia
# SatComm ex2
using DelimitedFiles
using Makie

const e = Base.MathConstants.e

function init()

    # ========== PartA ==========

    # ----- Input 1 -----
    percentage_of_time = [1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]
    regions_for_group4 = ['A', 'D', 'F', 'K']

    p837 = initP837()

    # ----- Input 2 -----
    link_availability = 99.99 # %
    index_of_interest = 1
    mind = 100
    for i=1:length(percentage_of_time)
        d = abs(100 - (percentage_of_time[i]+link_availability))
        if d < mind
            mind = d
            index_of_interest = i
        end
    end

    # ----- Input 3 -----
    frequency_range = collect(1:30) # GHz

    # ----- Input 4 -----
    elevation_angle1 = 30.0 # degrees
    elevation_angle2 = 80.0 # degrees

    # ----- Input 5 -----
    es_height = 0.6 # km

    # ----- Input 6 -----
    pol_tilt_angle = 30.0 # degrees

    φ = -60.0 # region of group4 (Yes, we are using a floating Earth station!)

    Ap_ang1 = calcAttenuation(percentage_of_time, regions_for_group4, p837, index_of_interest, frequency_range, elevation_angle1, pol_tilt_angle, es_height, φ, link_availability)

    Ap_ang2 = calcAttenuation(percentage_of_time, regions_for_group4, p837, index_of_interest, frequency_range, elevation_angle2, pol_tilt_angle, es_height, φ, link_availability)

    visualizeAttenuation(Ap_ang1, Ap_ang2, frequency_range)

    # ===========================

    # ========== PartB ==========

    frequency_range = collect(4:30)
    XPD = calcXPD(Ap_ang1, pol_tilt_angle, frequency_range, elevation_angle1, regions_for_group4, link_availability)

    visualizeXPD(XPD, frequency_range)

    # ===========================
end

function calcAttenuation(percentage_of_time::Array{Float64},
                regions::Array{Char},
                p837::Dict{Char, Array{Float64}},
                index_of_interest::Int64,
                frequency_range::Array{Int64},
                θ::Float64,
                τ::Float64,
                hs::Float64,
                φ::Float64,
                link_availability::Float64)

    Ap = Dict{Char, Array{Float64}}()
    Re = 8500 # km

    # Step 1
    h0 = initP839(φ)
    hr = h0 + 0.36 # km

    # Step 2
    for r_ in regions
        if (hr - hs) > 0
            Ls = 0
            if θ >= 5
                Ls = (hr - hs) / sind(θ)
            else
                Ls = 2 * (hr - hs) / (sqrt(sind(θ)^2 + (2 * (hr - hs) / Re)) + sind(θ))
            end

            # Step 3
            Lg = Ls * cosd(θ)

            # Step 4
            R = p837[r_][index_of_interest]
            if R == 0
                Ap[r_] = zeros(length(frequency_range))
                continue
            end
            # Step 5
            for f in frequency_range
                k, α = initP838(f, θ, τ)
                γR = k * R^α

                # Step 6
                r = 1 / (1 + 0.78 * sqrt(Lg * γR / f) - 0.38 * (1 - e^(-2 * Lg)))

                # Step 7
                Lr = 0
                ζ = atand((hr - hs) / (Lg * r))
                if ζ > θ
                    Lr = Lg * r / cosd(θ)
                else
                    Lr = (hr - hs) / sind(θ)
                end
                χ = 0
                if abs(φ) < 36
                    χ = 36 - abs(φ)
                end
                v = 1 / (1 + sqrt(sind(θ)) * (31 * (1 - e^(-θ / (1+χ)) * sqrt(Lr*γR) / f^2 - 0.45)))

                # Step 8
                Le = Lr * v

                # Step 9
                A = γR * Le

                # Step 10
                β = 0
                p = 100 - link_availability
                if p < 1 && abs(φ) < 36 && θ >= 25
                    β = -0.005 * (abs(φ) - 36)
                elseif p < 1 && abs(φ) < 36
                    β = -0.005 * (abs(φ) - 36) + 1.8 - 4.25 * sind(θ)
                end
                Ap_ = A * (p / 0.01)^(-(0.655 + 0.033 * log(p) - 0.045 * log(A) - β * (1 - p) * sind(θ)))
                if haskey(Ap, r_)
                    push!(Ap[r_], Ap_)
                else
                    Ap[r_] = [Ap_]
                end
            end
        else
            Ap[r_] = zeros(length(frequency_range))
        end
    end
    return Ap
end

function visualizeAttenuation(Ap1::Dict{Char, Array{Float64}}, Ap2::Dict{Char, Array{Float64}}, frequency_range::Array{Int64})
    println("Preparing graphics...")

    i = 1
    colours = [:yellow, :red, :green, :blue]
    regions = []

    # Dummy initialization of the two plots
    l = lines(zeros(length(frequency_range)), zeros(length(frequency_range)))
    l2 = lines(zeros(length(frequency_range)), zeros(length(frequency_range)))

    for a in Ap1
        region = a[1]
        ap = a[2]
        lines!(l, frequency_range, ap, color = colours[i])
        lines!(l2, frequency_range, Ap2[region], color = colours[i])
        push!(regions, a[1])
        i += 1
    end
    ax1 = l[Axis]
    ax1[:ticks][:title_gap] = 20
    ax1[:names][:axisnames] = ("frequency (GHz)", "30 degrees\nattenuation (dB)")
    ax2 = l2[Axis]
    ax2[:ticks][:title_gap] = 20
    ax2[:names][:axisnames] = ("frequency (GHz)", "80 degrees\nattenuation (dB)")
    #dummy initialization of the legend
    legend = text("_", color = :white, position = Vec(0,0), textsize = 0.0001, show_axis = false)
    for i=1:length(regions)
        text!(legend, "Region: "*regions[i], color = colours[i], position = Vec(0,i/2), textsize = 1, show_axis = false)
    end
    scene = vbox(hbox(l, l2), legend, sizes = [0.7, 0.3])
    display(scene)
    println("Press return to continue to the solution of Part B...")
    readline()
end

function calcXPD(Ap::Dict{Char, Array{Float64}}, τ::Float64, frequency_range::Array{Int64}, θ::Float64, regions::Array{Char}, link_availability::Float64)
    XPD = Dict{Char, Array{Float64}}()
    for region in regions
        for f in frequency_range
            # Step 1
            if f >= 6 && f <= 55
                Cf = 0
                if f >= 6 && f < 9
                    Cf = 60 * log10(f) - 28.3
                elseif f >= 9 && f < 36
                    Cf = 26 * log10(f) + 4.1
                else
                    Cf = 35.9 * log10(f) - 11.3
                end

                # Step 2
                Vf = 0
                if f >= 6 && f < 9
                    Vf = 30.8 * f^(-0.21)
                elseif f >= 9 && f < 20
                    Vf = 12.8 * f^0.19
                elseif f >= 20 && f < 40
                    Vf = 22.6
                elseif f >= 40 && f <= 55
                    Vf = 13.0 * f^0.15
                end
                # Here I am using the frequency as an index, which is generally 
                # dangerous, but in our current situation, we have f∈[4,30]
                Ca = Vf * log10(Ap[region][f])

                # Step 3
                Cτ = -10 * log10(1 - 0.484 * (1 + cosd(4 * τ)))

                # Step 4
                if θ <= 60
                    Cθ = -40 * log10(cosd(θ))
                # else what?!
                end

                # Step 5
                σ = initESD()
                p = round(100 - link_availability, sigdigits=4)
                Cσ = 0.0053 * σ[p]^2

                # Step 6
                XPDrain = Cf - Ca + Cτ + Cθ + Cσ

                # Step 7
                Cice = XPDrain * (0.3 + 0.1 * log10(p))/2

                # Step 8
                XPDp = XPDrain - Cice
                if haskey(XPD, region)
                    push!(XPD[region], XPDp)
                else
                    XPD[region] = [XPDp]
                end
            end
        end
        # Get frequencies smaller than 6
        f_ = reverse(frequency_range[frequency_range .< 6])
        validf = 
        for f in f_
            # Assuming at least one frequency higher than 6 was found, so that the XPD dict is not empty
            # Also assuming that the last element of the dict has the equivalent value for a maximum of 30GHz
            XPD1 = XPD[region][end]
            f1 = frequency_range[end]
            f2 = f
            XPD2 = XPD1 - 20 * log10((f2 * sqrt(1 - 0.484 * (1 + cosd(4 * τ)))) / (f1 * sqrt(1 - 0.484 * (1 + cosd(4 * τ)))))
            prepend!(XPD[region], XPD2)
        end
        println(XPD[region])
    end
    return XPD
end

function visualizeXPD(XPD::Dict{Char, Array{Float64}}, frequency_range::Array{Int64})
    println("Preparing graphics...")

    i = 1
    colours = [:yellow, :red, :green, :blue]
    regions = []

    # Dummy initialization of the two plots
    l = lines(zeros(length(frequency_range)), zeros(length(frequency_range)))

    for xpd_ in XPD
        region = xpd_[1]
        xpd = xpd_[2]
        lines!(l, frequency_range, xpd, color = colours[i])
        push!(regions, xpd_[1])
        i += 1
    end
    ax1 = l[Axis]
    ax1[:ticks][:title_gap] = 20
    ax1[:names][:axisnames] = ("frequency (GHz)", "30 degrees\n XPD (dB)")
    #dummy initialization of the legend
    legend = text("_", color = :white, position = Vec(0,0), textsize = 0.0001, show_axis = false)
    for i=1:length(regions)
        text!(legend, "Region: "*regions[i], color = colours[i], position = Vec(0,i/2), textsize = 1, show_axis = false)
    end
    scene = vbox(l, legend, sizes = [0.8, 0.2])
    display(scene)
    println("Press return to exit...")
    readline()
end


# Helper functions, nothing interesting here...
function initP837()
    p837 = Dict{Char, Array{Float64}}()
    p837['A'] = [0.1, 0.8, 2, 5, 8, 14, 22]
    p837['B'] = [0.5, 2, 3, 6, 12, 21, 32]
    p837['C'] = [0.7, 2.8, 5, 9, 15, 26, 42]
    p837['D'] = [2.1, 4.5, 8, 13, 19, 29, 42]
    p837['E'] = [0.6, 2.4, 6, 12, 22, 41, 70]
    p837['F'] = [1.7, 4.5, 8, 15, 28, 54, 78]
    p837['G'] = [3, 7, 12, 20, 30, 45, 65]
    p837['H'] = [2, 4, 10, 18, 32, 55, 83]
    p837['J'] = [8, 13, 20, 28, 35, 45, 55]
    p837['K'] = [1.5, 4.2, 12, 23, 42, 70, 100]
    p837['L'] = [2, 7, 15, 33, 60, 105, 150]
    p837['M'] = [4, 11, 22, 40, 63, 95, 120]
    p837['N'] = [5, 15, 35, 65, 95, 140, 180]
    p837['P'] = [12, 34, 65, 105, 145, 200, 250]
    p837['Q'] = [24, 49, 72, 96, 115, 142, 170]
    return p837
end

function initP838(frequency::Int64, θ::Float64, τ::Float64)
    kv = Dict{Int64, Float64}()
    αv = Dict{Int64, Float64}()
    kh = Dict{Int64, Float64}()
    αh = Dict{Int64, Float64}()

    kv[1] = 0.0000308
    αv[1] = 0.8592
    kv[2] = 0.0000998
    αv[2] = 0.9490
    kv[3] = 0.0001942
    αv[3] = 1.0688
    kv[4] = 0.0002461
    αv[4] = 1.2476
    kv[5] = 0.0002428
    αv[5] = 1.5317
    kv[6] = 0.0004878
    αv[6] = 1.5728
    kv[7] = 0.001425
    αv[7] = 1.4745
    kv[8] = 0.003450
    αv[8] = 1.3797
    kv[9] = 0.006691
    αv[9] = 1.2895
    kv[10] = 0.01129
    αv[10] = 1.2156
    kv[11] = 0.01731
    αv[11] = 1.1617
    kv[12] = 0.02455
    αv[12] = 1.1216
    kv[13] = 0.03266
    αv[13] = 1.0901
    kv[14] = 0.04126
    αv[14] = 1.0646
    kv[15] = 0.05008
    αv[15] = 1.0440
    kv[16] = 0.05899
    αv[16] = 1.0273
    kv[17] = 0.06797
    αv[17] = 1.0137
    kv[18] = 0.07708
    αv[18] = 1.0025
    kv[19] = 0.08642
    αv[19] = 0.9930
    kv[20] = 0.09611
    αv[20] = 0.9847
    kv[21] = 0.1063
    αv[21] = 0.9771
    kv[22] = 0.1170
    αv[22] = 0.9700
    kv[23] = 0.1284
    αv[23] = 0.9630
    kv[24] = 0.1404
    αv[24] = 0.9561
    kv[25] = 0.1533
    αv[25] = 0.9491
    kv[26] = 0.1669
    αv[26] = 0.9421
    kv[27] = 0.1813
    αv[27] = 0.9349
    kv[28] = 0.1964
    αv[28] = 0.9277
    kv[29] = 0.2124
    αv[29] = 0.9203
    αv[30] = 0.2291
    kv[30] = 0.9129

    kh[1] = 0.0000259
    αh[1] = 0.9691
    kh[2] = 0.0000847
    αh[2] = 1.0664
    kh[3] = 0.0001390
    αh[3] = 1.2322
    kh[4] = 0.0001071
    αh[4] = 1.6009
    kh[5] = 0.0002162
    αh[5] = 1.6969
    kh[6] = 0.0007056
    αh[6] = 1.5900
    kh[7] = 0.001915
    αh[7] = 1.4810
    kh[8] = 0.004115
    αh[8] = 1.3905
    kh[9] = 0.007535
    αh[9] = 1.3155
    kh[10] = 0.01217
    αh[10] = 1.2571
    kh[11] = 0.01772
    αh[11] = 1.2140
    kh[12] = 0.02386
    αh[12] = 1.1825
    kh[13] = 0.03041
    αh[13] = 1.1586
    kh[14] = 0.03738
    αh[14] = 1.1396
    kh[15] = 0.04481
    αh[15] = 1.1233
    kh[16] = 0.05282
    αh[16] = 1.1086
    kh[17] = 0.06146
    αh[17] = 1.0949
    kh[18] = 0.07078
    αh[18] = 1.0818
    kh[19] = 0.08084
    αh[19] = 1.0691
    kh[20] = 0.09164
    αh[20] = 1.0568
    kh[21] = 0.1032
    αh[21] = 1.0447
    kh[22] = 0.1155
    αh[22] = 1.0329
    kh[23] = 0.1286
    αh[23] = 1.0214
    kh[24] = 0.1425
    αh[24] = 1.0101
    kh[25] = 0.1571
    αh[25] = 0.9991
    kh[26] = 0.1724
    αh[26] = 0.9884
    kh[27] = 0.1884
    αh[27] = 0.9780
    kh[28] = 0.2051
    αh[28] = 0.9679
    kh[29] = 0.2224
    αh[29] = 0.9580
    αh[30] = 0.2403
    kh[30] = 0.9485

    k = (kh[frequency] + kv[frequency] + (kh[frequency] - kv[frequency]) * cosd(θ)^2 * cosd(2 * τ)) / 2
    α = (kh[frequency] * αh[frequency] + kv[frequency] * αv[frequency] + (kh[frequency] * αh[frequency] - kv[frequency] * αv[frequency]) * cosd(θ)^2 * cosd(2 * τ)) / (2 * k)

    return k, α
end

# Effective Standard Deviation
function initESD()
    σ = Dict{Float64, Int64}()
    σ[1.0] = 0
    σ[0.1] = 5
    σ[0.01] = 10
    σ[0.001] = 15
    return σ
end

# Isotherm height
function initP839(lat::Float64)
    lats = readdlm(String(@__DIR__)*"/ITU-R/Lat.txt")
    h0s = readdlm(String(@__DIR__)*"/ITU-R/h0.txt")
    for i=1:size(lats)[1]
        for j=1:size(lats)[2]
            if lats[i,j] == lat
                return h0s[i,j]
            end
        end
    end
    return -1
end

init()
