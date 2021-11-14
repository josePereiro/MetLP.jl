struct FBAOut{T<:Real} <: AbstractMetState
    v::Vector{T} # Flux vector solution
    obj_val::T # The value of the objective
    obj_ider # The used objective identifier
    sol # The LP solution object 
end

# empty FBAOut
objval(out::FBAOut) = out.obj_val
objider(out::FBAOut) = out.obj_ider
LPsol(out::FBAOut) = out.sol