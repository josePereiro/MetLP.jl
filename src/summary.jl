import MetNets._print_state_head

function _print_state_head(io::IO, out::FBAOut)
    printstyled(io, " FBAOut: objective val: $(out.obj_val)\n", color = MetNets.INFO_COLOR)
end
