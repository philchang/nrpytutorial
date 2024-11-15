"""
Some simple scripts needed for MaNGa

Author: Thiago Assumpção; assumpcaothiago@gmail.com
"""

import nrpy.indexedexp as ixp
import nrpy.helpers.jacobians as jac
import nrpy.c_codegen as ccg

T4CartUU = ixp.declarerank2("T4CartUU", dimension=4, symmetry="01")

CoordSystem = "Spherical"

T4UU = jac.basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis(
    CoordSystem, T4CartUU
)

body = ccg.c_codegen(
    [T4UU[i][j] for i in range(4) for j in range(i, 4)],
    [f"T4UU{i}{j}" for i in range(4) for j in range(i, 4)],
    verbose=False,
    include_braces=False,
)

if __name__ == "__main__":
    print("Generating file: transform_T4UU_from_Cart_to_spherical.h")
    with open("transform_T4UU_from_Cart_to_spherical.h", "w") as file:
        file.write(body)
