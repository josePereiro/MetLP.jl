import Base.show
import Base.isnan
import Base.isempty

show(io::IO, out::FBAOut) = _print_state_head(io, out)

isnan(out::FBAOut)  = isnan(objval(out))
isempty(out::FBAOut)  = isnan(out)