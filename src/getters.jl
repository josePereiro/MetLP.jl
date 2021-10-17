# MetState interface
av(out::FBAOut) = out.v
va(out::FBAOut) = zeros(length(out.v))