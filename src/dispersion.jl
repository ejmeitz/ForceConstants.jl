export get_dispersion_points

function get_dispersion_points(sys::UnitCellSystem{3}, pot::PairPotential; directions = ([1.0, 0.0, 0.0],);
                                    group_by = :branch, tol = 1e-12, unit_system = :REAL)
    
    @assert group_by ∈ [:branch, :k_point] "group_by can only be :branch or :k_point"

    Δk = (2*pi)./sys.box_sizes_SC
    left_zone_edges = (-pi./sys.box_sizes_UC) .+ Δk
    right_zone_edges = pi./sys.box_sizes_UC

    k_points = []
    k_current = ustrip.([left_zone_edges[1],0,0])
        # direction ./= norm(direction) #convert to unit vector
        
    #This is just 1,0,0 direction for now
    for i in range(left_zone_edges[1], right_zone_edges[1], step = Δk[1])
        push!(k_points, copy(k_current))
        k_current .+= [ustrip(Δk[1]),0,0]
    end

    ω_all = Dict()
    println(k_points)
    for k_point in k_points
        dynmat_uc = dynamicalMatrix(sys, pot, SVector{3}(k_point), tol)
        ω_sq = eigen(Hermitian(ustrip.(dynmat_uc))).values;
        ω_all[k_point] = ω_sq
    end
    #Sqrt and convert units on freqs
    for k in keys(ω_all)

        if ustrip.(k) == [0.0, 0.0, 0.0] #if at origin set rigid translation modes to 0
            idxs = sortperm(abs.(ω_all[k]))
            ω_all[k][idxs[1:3]] .= 0.0
        end

        if unit_system == :REAL
            ω_all[k] = sqrt.(ω_all[k] .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)
        else
            throw(error("Only real units supported at this time"))
        end
    end

    return ω_all
end