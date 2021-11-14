import Base.isnan
import Base.isempty

isnan(out::FBAOut)  = isnan(objval(out))
isempty(out::FBAOut) = isnan(out)