import Base.show

show(io::IO, out::FBAOut) = _print_state_head(io, out)