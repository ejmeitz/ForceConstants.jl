
#could remove allocations here and just make it a ! function
function nearest_mirror(r_ij, box_sizes)
    r_x = r_ij[1]; r_y = r_ij[2]; r_z = r_ij[3]
    L_x, L_y, L_z = box_sizes
    if r_x > L_x/2
        r_x -= L_x
    elseif r_x < -L_x/2
        r_x += L_x
    end
        
    if r_y > L_y/2
        r_y -= L_y
    elseif r_y < -L_y/2
        r_y += L_y  
    end
        
    if r_z > L_z/2
        r_z -= L_z
    elseif r_z < -L_z/2
        r_z += L_z
    end
    
    return [r_x,r_y,r_z] 
end