module Select
    include("select/_select.jl")
    
    export values_along, values_in, indices_and_weights, patches_along, 
    patch_at, is_inside, corner_indices, patches_in, box, plane, patches_in, corner_indices_all, values_along2
end