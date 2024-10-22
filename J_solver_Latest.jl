module J_solver

include( "./ConstantsAndFunctions.jl" )
using .ConstantsAndFunctions
using ProgressMeter

# Generate r, nHI, and T arrays
Nr    = 1000
r     = LinRange(0, 1*pc, Nr)      # r from 0 to 1 pc
nHI   = fill(1e2, Nr)              # Constant nHI
T     = fill(1e2, Nr)              # Constant T


rsbar = 0.05
L_Lyα = 1.0*Lsun

#println( nHI .* σ0.(T) )

""" Generates a suitable frequency grid for a given set of
    cloud properties.  """
function GetFrequencies(r::AbstractVector{Float64}, nHI::AbstractVector{Float64}, 
                              T::AbstractVector{Float64}, Ncoarse = 400, ηfreq = 4.0)
    # Compute τ_cl:
    τcl = τ(r, nHI, T)

    # Get maximum and minimum temperature, and corresponding a(T)
    # for the latter:
    Tmax = maximum( T )
    Tmin = minimum( T ); aTmin = a(Tmin)

    # Get the minimum and maximum Doppler width:
    Δν_Dmin = (b(Tmin)/c)*νLyα
    Δν_Dmax = (b(Tmax)/c)*νLyα

    # Generate frequency grid that resolves down < Δν_Dmin around
    # line center:
    ν_finegrid_min   = νLyα - 10.0*Δν_Dmax                   # Maximum frequency with fine gridding
    ν_finegrid_max   = νLyα + 10.0*Δν_Dmax                   # Minimum frequency with fine gridding
    Δν_fine          = 0.1*Δν_Dmin                           # Spacing of fine gridding
    ν_finegrid       = ν_finegrid_min:Δν_fine:ν_finegrid_max # Fine grid near line center

    # Generate coarser gridding outside:
    Δν_max           = ηfreq * Δν_Dmax * ( aTmin * τcl )^(1/3) # Maximum considered frequency width
    ν_coarsegrid_min = ν_finegrid_min - Δν_max
    ν_coarsegrid_max = ν_finegrid_max + Δν_max

    ν_coarsegrid_low = LinRange( ν_coarsegrid_min, ν_finegrid_min - Δν_fine, div( Ncoarse, 2 ) )
    ν_coarsegrid_up  = LinRange( ν_finegrid_max + Δν_fine, ν_coarsegrid_max, div( Ncoarse, 2 ) )

    # Put together to get the full frequency grid:
    ν = collect( vcat( ν_coarsegrid_low, ν_finegrid, ν_coarsegrid_up ) )

    return ν
end





""" Calculates all the needed factors that depend on cloud conditions,
    needed to numerically solve for J(r,ν). 
        
    Arguments:
    
    r     :   Array of radial grid points (cm)
    ν     :   Array of frequency grid points (Hz)
    nHI   :   Array of HI number densities at radial grid points (cm^-3)
    T     :   Array of gas temperatures at radial grid points (K)
    rsbar :   Ratio of emission region to cloud radius (minimum value set automatically by r[1])
    L_Lyα :   Total Lyα luminosity from sources (erg/s) 
    
    Returns following arrays: 
    
    R_{i+1/2}, R_{i-1/2}, F_{i,j+1/2}, F_{i,j-1/2}
    
    """
function GetFactors(r::AbstractVector{Float64}, ν::AbstractVector{Float64}, nHI::AbstractVector{Float64}, 
                    T::AbstractVector{Float64}, rsbar::Float64, L_Lyα::Float64)
    Nr = length(r)
    Nν = length(ν)

    # Initialize arrays:

    Ri_pj   = zeros(Nr, Nν)
    Ri_mj   = zeros(Nr, Nν)
    Fi_jp   = zeros(Nr, Nν)
    Fi_jm   = zeros(Nr, Nν)

    # Calculate midpoints of r and ν:
    rmid = (r[1:Nr-1] .+ r[2:Nr])./2 # Midpoint values r_{i±1/2}  
    νmid = (ν[1:Nν-1] .+ ν[2:Nν])./2 # Midpoint values ν_{j±1/2}

    Δr = zeros(Nr)
    for i in 2:Nr-1
        Δr[i] = (r[i+1] - r[i-1])/2
    end
    Δν = zeros(Nν)
    for j in 2:Nν-1
        Δν[j] = (ν[j+1] - ν[j-1])/2
    end

    # Calculate α = n_HI*σ0 at each radial point & mid-point:
    α    = nHI .* σ0.(T)
    αmid = (α[1:Nr-1] .+ α[2:Nr])./2 # Midpoint values α_{i±1/2}

    # Calculate Δν_D at every radial point:
    Δν_D = (b.(T) ./ c) .* νLyα

    # Calculate H(x,a) at every grid point and interface:

    Hij  = zeros(Nr, Nν)
    Hijp = zeros(Nr, Nν)
    Hijm = zeros(Nr, Nν)
    for i in 1:Nr
        for j = 1:Nν
            #Hij[i,j]  = H( (ν[j] - νLyα)/Δν_D[i], a(T[i]) )
            #Hijp[i,j] = H( (ν[j] + Δν[j]/2 - νLyα)/Δν_D[i], a(T[i]) ) # H_{i,j+1/2}
            #Hijm[i,j] = H( (ν[j] - Δν[j]/2 - νLyα)/Δν_D[i], a(T[i]) ) # H_{i,j-1/2}
            Δx = Δν[j]/Δν_D[i]
            Hij[i,j] = Htophat( (ν[j] - νLyα)/Δν_D[i],  a(T[i]),  Δx )
        end
    end

    ### GET RADIAL FLUX PRE-FACTORS ###
    for i = 2:Nr-1
        for j = 2:Nν-1
            Himidpj   = (Hij[i,j] + Hij[i+1,j])/2  # H_{i+1/2,j}
            Himidmj   = (Hij[i,j] + Hij[i-1,j])/2  # H_{i-1/2,j}

            # Volume element ΔV_i = (4π/3)*(r_{i+1/2}^3 - r_{i-1/2}^3) :
            ΔV        = (4*π/3)*( rmid[i]^3 - rmid[i-1]^3 )
            
            # Radial flux pre-factors:
            Ri_pj[i,j] = 8 * π * rmid[i]^2 / (ΔV * αmid[i] * Himidpj * (Δr[i] + Δr[i+1]))      # R_{i+1/2,j}
            Ri_mj[i,j] = 8 * π * rmid[i-1]^2 / (ΔV * αmid[i-1] * Himidmj * (Δr[i] + Δr[i-1]))  # R_{i-1/2,j}
        end
    end

    ### GET FREQUENCY DIFFUSION PRE-FACTORS & EMISSION FACTORS ###

    # Get nearest radius corresponding to rsbar, and source region volume:

    Rs_exact     = rsbar * maximum(r)
    source_index = max( 2, argmin( abs.(r .- Rs_exact) ) )
    Rs           = r[source_index]  # Closest discrete source radius
    Vs           = (4*π/3)*(Rs^3)   # Source region volume

    # Lyα emissivity:

    js  = zeros( Nr, Nν )

    Tem    = 1e3
    Δν_Dem = (b(Tem)/c)*νLyα

    xem  = (ν .- νLyα)./Δν_Dem
    ϕ_ν = H.( xem, a( Tem ) )./( sqrt(π) .* Δν_Dem)

    for i in 1:Nr-1
        for j in 2:Nν-1
            #js[i,j] = L_Lyα*ϕ_ν[j] /(4*π*Vs)
            n = 10.0
            F = exp( - (r[i]/r[source_index])^n )
            js[i,j] = L_Lyα*Hij[i,j]*F/(4*π*Vs)
        end
    end

    for i = 1:Nr
        for j = 2:Nν-1
            #f = 0.65
            #Hijp = f*Hij[i,j] + (1-f)*Hij[i,j+1]
            #Hijp = min( Hij[i,j], Hij[i,j+1] )
            #Hijm = f*Hij[i,j] + (1-f)*Hij[i,j-1]
            #Hijm = min( Hij[i,j], Hij[i,j-1] )

            # Frequency diffusion factors F:
            Fi_jp[i, j] =  3 * α[i] * (Δν_D[i]^2) / (Δν[j]*(Δν[j]/Hij[i,j] + Δν[j+1]/Hij[i,j+1]))  # F_{i,j+1/2}
            Fi_jm[i, j] =  3 * α[i] * (Δν_D[i]^2) / (Δν[j]*(Δν[j]/Hij[i,j] + Δν[j-1]/Hij[i,j-1]))  # F_{i,j-1/2}

        end
    end


    return Ri_pj, Ri_mj, Fi_jp, Fi_jm, js, αmid, Hij
end

#ν =  GetFrequencies(r, nHI, T)
#Ri_pj, Ri_mj, Fi_jp, Fi_jm, js, αmid, Hmidedge = GetFactors(r, ν, nHI, T, rsbar, L_Lyα) 

#println( Fi_jm[2,:2] )


""" The solver for J(r,ν). Uses an iterative relaxation method 
    (SOR or Gauss-Seidel, depending on ω). Arguments:

    r     :   Array of radial grid points (cm)
    ν     :   Array of frequency grid points (Hz)
    nHI   :   Array of HI number densities at radial grid points (cm^-3)
    T     :   Array of gas temperatures at radial grid points (K)
    rsbar :   Ratio of emission region to cloud radius (minimum value set automatically by r[1])
    L_Lyα :   Total Lyα luminosity from sources (erg/s)
    
    Returns:
    
    J """
function J_sol(r::AbstractVector{Float64}, nHI::AbstractVector{Float64}, T::AbstractVector{Float64}, 
               rsbar::Float64, L_Lyα::Float64, spectrum = true, 
               max_iter = 20000, rtol = 1e-3, ω = 1.93)

    # Generate frequency grid for the given cloud:
    ν  = GetFrequencies( r, nHI, T )

    # Get grid numbers:
    Nr = length(r); Nν = length(ν)

    # Compute all needed pre-factors:
    Ri_pj, Ri_mj, Fi_jp, Fi_jm, js, αmid, Hij = GetFactors(r, ν, nHI, T, rsbar, L_Lyα)

    # Initialize J_ij ≡ J(r_i, ν_j) & error matrix:
    J      =   zeros( Nr, Nν ) 
    error  =   ones( Nr, Nν )
    error[:,1] .= 0.0; error[:,Nν] .= 0.0
    error[1,:] .= 0.0; error[Nr,:] .= 0.0

    #### SOR/GAUSS-SEIDEL ITERATION ####

    println("Starting iteration...")
    println("")
    @showprogress 1 "Solving for J:" for iter = 1:max_iter
        for i in 2:Nr-1
            for j in 2:Nν-1
                # Standard Gauss-Seidel update the intensity J[i, j]
                J_new    = (3*js[i,j] + Ri_pj[i,j]*J[i+1,j] + Ri_mj[i,j]*J[i-1,j] + Fi_jp[i,j]*J[i,j+1] + Fi_jm[i,j]*J[i,j-1]) /
                           (Ri_pj[i,j] + Ri_mj[i,j] + Fi_jp[i,j] + Fi_jm[i,j])

                # Use successive over-relaxation (if ω ≠ 1):
                J_new = ω*J_new + (1 - ω)*J[i,j]
                        
                # Calculate the difference for convergence checking
                error[i,j] = abs(J_new - J[i,j])/abs(J[i,j])

                # Update J to J_new:
                J[i, j]  = J_new 
            end
        end

        #### ENFORCE BOUNDARY CONDITIONS ######

        for j in 2:Nν-1
            # No flux at cloud center:
            J[1, j]   = J[2, j]

            #R1 = 6 / (αmid[1] * (r[2]-r[1])^2)

            #J[1, j] = ( 3*Hjs[1,j] + R1*J[2,j] + Fi_jp[1,j]*J[1,j+1] + Fi_jm[1,j]*J[1,j-1] ) /
            #          ( R1 + Fi_jp[1,j] + Fi_jm[1,j] )
            
            # No incoming intensity at cloud edge:
            #Δr_edge   = r[Nr] - r[Nr-1]
            #J[end, j] = J[end-1, j] / (1 + sqrt(3)*αmid[end-1]*Hmidedge[j]*Δr_edge)
            #J[Nr, j]  = J[Nr-1,j] / (1 + sqrt(3)*αmid[Nr-1]*Hij[Nr,j]*Δr_edge)

            #RN_p     = (r[Nr] + Δr_edge)^2 / ( r[Nr]^2 * Δr_edge^2 * αmid[Nr-1] )
            #Δτedge   = 2 * sqrt(3) * αmid[Nr-1] * Δr_edge

            #J[Nr, j] = ( 3*Hjs[Nr,j] + (RN_p + Ri_m[Nr])*J[Nr-1,j] + Fi_jp[Nr,j]*J[Nr,j+1] + Fi_jm[Nr,j]*J[Nr,j-1]) /
            #           ((1+Δτedge)*RN_p + Ri_m[Nr-1] + Fi_jp[Nr,j] + Fi_jm[Nr,j])
            J[Nr, j] = 0.0
        end

        # No intensity far from line center:
        J[:, 1]  .= 0.0 # J(ν << νLyα) = 0
        J[:, Nν] .= 0.0 # J(ν >> νLyα) = 0

        #### CHECK CONVERGENCE ###############

        max_error = maximum( error )

        if iter > min( 300, max_iter-1 ) && max_error < rtol
            println("")
            println("Converged successfully after $iter iterations!")
            break
        elseif iter == max_iter && max_error >= rtol
            println("")
            println("WARNING: Did not converge within error tolerance.")
        end

    end

    if spectrum == true
        # Compute emergent normalized spectrum
        Jedge = J[Nr-1, :]             # J at the cloud edge
        Δν_D  = (b(T[Nr])/c)*νLyα    # Doppler width at cloud edge
        x     = (ν .- νLyα)/Δν_D     # Dimensionless frequency

        # Normalize Jedge:
        Δx             = diff(x)
        Integral_Jedge = sum( Δx .* (Jedge[1:Nν-1] .+ Jedge[2:Nν]) ) / 2.0

        # Normalize Jedge:
        Jedge_normalized = Jedge / Integral_Jedge 
    end

    Δν = diff(ν); Δr = diff(r)
    volumes = 4 .* π .*(r[1:Nr-1] .^2) .* Δr

    # Frequency integral ∫dν J(r,ν) for each radial point
    Jtot = zeros(Nr)
    for i in 1:Nr
        Jtot[i] = sum(0.5 * (J[i, 1:Nν-1] .+ J[i, 2:Nν]) .* Δν)
    end
    
    # Total Lyα energy within the cloud:
    #E_Lyα = (4*π/c) * sum(Jtot .* volumes)

    # Get trapping time:
    #t_trap = (E_Lyα/L_Lyα)/( r[Nr]/c )

    return J, x, Jedge_normalized, Jtot
    #return  J, t_trap, Jtot
end

using PyPlot

# Assume J_sol is your current function and you have already executed it to get the results.
J, x, Jedge_normalized, Jtot = J_sol(r, nHI, T, rsbar, L_Lyα)

# Normalize r by r_max
r_norm = r ./ maximum(r)

# Create the plot
fig, ax = subplots()

# Use imshow to plot J as a function of x and r/rmax
cplot = ax.imshow(log10.(J), extent=(minimum(x), maximum(x), minimum(r_norm), maximum(r_norm)), 
              aspect="auto", origin="lower", cmap="viridis")

# Add colorbar
fig.colorbar(cplot, ax=ax)

# Labels and title
xlabel("x = (ν - ν_Lyα) / Δν_D")
ylabel("r/R_{cl}")
title("Intensity J as a function of x and r / r_max")

#plot( x, Jedge_normalized )

# Show the plot
show()



#J_sol()
end # End of module