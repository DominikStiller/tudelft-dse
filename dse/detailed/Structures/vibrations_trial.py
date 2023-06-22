from material_properties import materials

# Let's try to do a 2m long, square beam
material = materials["Al/Si"]
L = 2
side = 5e-2
A = side**2
I = side**4 / 12
E = material.E
mass = A * L * material.rho
