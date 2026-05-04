using BlockArrays, SparseArrays
order = 2
larger_order = 2order
basis_length = larger_order+1

in_order = larger_order
out_order = order

U_arr = [ 
    [dl dm 0.0im; l m 0.0im; 0.0im 0.0im 0.0im]
for dl = 0:larger_order for dm = -dl:dl for l = 0:order for m = -l:l];

Ublocks = reshape(U_arr, ((order+1)^2, (larger_order+1)^2))
V = sparse(mortar(Ublocks))

V[end-2:end,1:3]


order = 1
ω = 1.0;
3 * (order + 1)^2
x = [1.0,1.1,1.2]
vs_in = regular_basis_function(medium, ω, field_type)(order,x)
us_in = outgoing_basis_function(medium, ω, field_type)(order,x)
us_in[:,1:3]

us_in[:,4:6]
us_in[:,7:9]
us_in[:,10:12]

vs_in[:,4:6]

order = 1
pbasis = pressure_regular_basis(ω, x, medium, order, field_type) 
Φbasis = shearΦ_regular_basis(ω, x, medium, order, field_type) 
χbasis = shearχ_regular_basis(ω, x, medium, order, field_type) 


reshape([pbasis Φbasis χbasis] |> transpose,3,:) - ([pbasis; Φbasis; χbasis] |> transpose)