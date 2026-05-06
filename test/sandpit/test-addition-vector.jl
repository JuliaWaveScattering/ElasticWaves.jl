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

medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0)
order = 1
ω = 1.0;
3 * (order + 1)^2
x = [1.0,1.1,1.2]
vs_in = regular_basis_function(medium, ω, field_type)(order,x)
vs_in[:,4:6]

order = 1
pbasis = pressure_regular_basis(ω, x, medium, order, field_type) 
Φbasis = shearΦ_regular_basis(ω, x, medium, order, field_type) 
χbasis = shearχ_regular_basis(ω, x, medium, order, field_type) 

reshape([pbasis Φbasis χbasis] |> transpose,3,:) - vs_in |> norm

ω = 1.0;
r = [-0.16760003217875297,-0.35146260115553596,-0.2124509700641245];
d = [-2.9010760893418666, -3.353041498826335,0.12643092521032334];

order = 1;
vs_in = regular_basis_function(medium, ω, field_type)(order,r)




V = regular_translation_matrix(medium, 1, 1, ω, d, field_type)

V*transpose(vs_in)

V[1:3,:] |> collect 

V[1:3,4:end] |> collect
V[4:6,4:end] |> collect 