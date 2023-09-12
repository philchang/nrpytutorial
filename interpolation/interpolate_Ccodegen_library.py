from outputC import add_to_Cfunction_dict

def indent_Ccode(indent, Ccode):
    Ccodesplit = Ccode.splitlines()
    outstring = ""
    for i in range(len(Ccodesplit)):
        outstring += indent + Ccodesplit[i] + '\n'
    return outstring

def desc_format(desc):
    return "/*\n" + indent_Ccode(" * ", desc) + " */\n"

# We just return a string containing the function so the function can be more easily inlined.
#  Otherwise (just_return_string=False) will result in function being placed in its own C file.
def codegen_Lagrange_interp_coeffs_3D(just_return_string=True):
    includes = ["NRPy_basic_defines.h", # // Defines REAL, NinterpGHOSTS, NGHOSTS, & interpstruct. This header is generated in NRPy+.
                "stdio.h", "stdlib.h", "math.h"]
    desc = "Sets up interpolation coefficients at a given point"
    name = "Lagrange_interp_coeffs_3D"
    params = """const int N0,const int N1,const int N2,
                               const REAL *restrict x0i,const REAL *restrict x1j,const REAL *restrict x2k,
                               const REAL x012[3],
                               REAL *restrict l0i__times__w0i_inv,REAL *restrict m1j__times__w1j_inv,REAL *restrict n2k__times__w2k_inv"""
    body = r"""
// Computing the denominator using integer operations in this way
//     yields 6% overall speedup, with 128^3 grid,
//     num_interp_pts=1e7, and poly interp order=4
// Total cost: [INTEGER: (N0 + N1 + N2) * (1 assign + 1 mul + 1 sub) + 6 assigns] + 3 typecasts
for(int i=0;i<=N0;i++) {
    REAL prod_numer_i = 1.0; // l0i
    int  prod_denom_i = 1;   // w0i
    for(int l=0;  l<i;  l++) { prod_denom_i *= i-l; prod_numer_i *= x012[0] - x0i[l]; }
    for(int l=i+1;l<=N0;l++) { prod_denom_i *= i-l; prod_numer_i *= x012[0] - x0i[l]; }
    //                         [       w0i        ] [            l0i                ]
    l0i__times__w0i_inv[i] = prod_numer_i / ( (REAL)prod_denom_i );
}
for(int j=0;j<=N1;j++) {
    REAL prod_numer_j = 1.0; // m1j
    int  prod_denom_j = 1;   // w1j
    for(int m=0;  m<j;  m++) { prod_denom_j *= j-m; prod_numer_j *= x012[1] - x1j[m]; }
    for(int m=j+1;m<=N1;m++) { prod_denom_j *= j-m; prod_numer_j *= x012[1] - x1j[m]; }
    //                         [       w1j        ] [            m1j                ]
    m1j__times__w1j_inv[j] = prod_numer_j / ( (REAL)prod_denom_j );
}
for(int k=0;k<=N2;k++) {
    REAL prod_numer_k = 1.0; // n2k
    int  prod_denom_k = 1;   // w2k
    for(int n=0;  n<k;  n++) { prod_denom_k *= k-n; prod_numer_k *= x012[2] - x2k[n]; }
    for(int n=k+1;n<=N2;n++) { prod_denom_k *= k-n; prod_numer_k *= x012[2] - x2k[n]; }
    //                         [       w2k        ] [            n2k                ]
    n2k__times__w2k_inv[k] = prod_numer_k / ( (REAL)prod_denom_k );
}
"""
    body = indent_Ccode('    ', body)
    if just_return_string:
        return desc_format(desc) + "void" + " " + name + "(" + params + ") {\n" + body + "}\n"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body, enableCparameters=False)

# We just return a string containing the function so the function can be more easily inlined.
#  Otherwise (just_return_string=False) will result in function being placed in its own C file.
def codegen_Lagrange_sum_3D(just_return_string=True):
    includes = ["NRPy_basic_defines.h", # // Defines REAL, NinterpGHOSTS, NGHOSTS, & interpstruct. This header is generated in NRPy+.
                "stdio.h", "stdlib.h", "math.h"]
    desc = """Computes 3D sum for interpolation at a single point.
    Note that this function does NOT multiply by the requisite (Delta x_i)^(N_i)"""
    type = "inline REAL"
    name = "Lagrange_sum_3D"
    params = """const int N0,const int N1,const int N2,
    const REAL *restrict l0i__times__w0i_inv,
    const REAL *restrict m1j__times__w1j_inv,
    const REAL *restrict n2k__times__w2k_inv,
    const REAL *restrict f"""
    body = r"""
// Now perform sum, with a total cost of
//  N0N1 * (1 assignment, 1 multiply [1 + 1 ~ 2 FLOPs])
//   + N0N1N2 * ( 1 assignment, 1 add, 2 multiplies [1 + 1 + 2 ~ 4 FLOPs])
REAL sum = 0.0;
int idx = 0;
for(int i=0;i<=N0;i++) {
    for(int j=0;j<=N1;j++) {
        // (N0+1)(N1+1) * (1 assignment, 1 multiply):
        const REAL l0i_term_times_m1j_term = l0i__times__w0i_inv[i]*m1j__times__w1j_inv[j];
        for(int k=0;k<=N2;k++) {
            // (N0+1)(N1+1)(N2+1) * ( 1 assignment, 1 add, 2 multiplies ):
            sum += f[idx] * ( l0i_term_times_m1j_term *  n2k__times__w2k_inv[k] );
            idx++;
        }
    }
}
return sum;
"""
    body = indent_Ccode('    ', body)
    if just_return_string:
        return desc_format(desc) + type + " " + name + "(" + params + ") {\n" + body + "}\n"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        type=type,     name=name,     params=params,
        body=body, enableCparameters=False)

# We just return a string containing the function so the function can be more easily inlined.
#  Otherwise (just_return_string=False) will result in function being placed in its own C file.
def codegen_construct_x0i_x1j_x2k_stencil_arrays_3D(just_return_string=True):
    includes = ["NRPy_basic_defines.h", # // Defines REAL, NinterpGHOSTS, NGHOSTS, & interpstruct. This header is generated in NRPy+.
                "stdio.h", "stdlib.h", "math.h"]
    desc = """Construct integer arrays x0i_idx[N0],x1j_idx[N1],x2k_idx[N2] """
    name = "construct_x0i_x1j_x2k_stencil_arrays_3D"
    params = """const int Nxx_plus_2NGHOSTS0,const int Nxx_plus_2NGHOSTS1,const int Nxx_plus_2NGHOSTS2,  REAL *restrict xx[3],
    const REAL x012[3],
    const int N0,const int N1,const int N2,
    int *restrict x0i_idx,int *restrict x1j_idx,int *restrict x2k_idx"""
    body = r"""
// Total cost: 3 * (2 subs + 1 div + 1 add) + 3 * (1 integer sub, div, and add)
    const int i012_int[3] = { (int)( (x012[0] - xx[0][0])/(xx[0][1] - xx[0][0]) + 0.5 ),
                              (int)( (x012[1] - xx[1][0])/(xx[1][1] - xx[1][0]) + 0.5 ),
                              (int)( (x012[2] - xx[2][0])/(xx[2][1] - xx[2][0]) + 0.5 ) };
//#define BOUNDS_CHECKING
#ifndef BOUNDS_CHECKING
    for(int i0=0;i0<=N0;i0++) x0i_idx[i0] = i012_int[0] - N0/2 + i0;
    for(int i1=0;i1<=N1;i1++) x1j_idx[i1] = i012_int[1] - N1/2 + i1;
    for(int i2=0;i2<=N2;i2++) x2k_idx[i2] = i012_int[2] - N2/2 + i2;
#else
    //printf("%d %d %d | %d iiii\n", i012_int[0], i012_int[1], i012_int[2], Nxx_plus_2NGHOSTS0);
    for(int i0=0;i0<=N0;i0++) {
        x0i_idx[i0] = i012_int[0] - N0/2 + i0;
        if(x0i_idx[i0] < 0 || x0i_idx[i0] >= Nxx_plus_2NGHOSTS0) {
            printf("ERROR: index %d is out of range!\n",x0i_idx[i0]); exit(1);
        }
    }
    for(int i1=0;i1<=N1;i1++) {
        x1j_idx[i1] = i012_int[1] - N1/2 + i1;
        if(x1j_idx[i1] < 0 || x1j_idx[i1] >= Nxx_plus_2NGHOSTS1) {
            printf("ERROR: index %d is out of range!\n",x1j_idx[i1]); exit(1);
        }
    }
    for(int i2=0;i2<=N2;i2++) {
        x2k_idx[i2] = i012_int[2] - N2/2 + i2;
        if(x2k_idx[i2] < 0 || x2k_idx[i2] >= Nxx_plus_2NGHOSTS2) {
            printf("ERROR: index %d is out of range!\n",x2k_idx[i2]); exit(1);
        }
    }
#endif
"""
    body = indent_Ccode('    ', body)
    if just_return_string:
        return desc_format(desc) + "void" + " " + name + "(" + params + ") {\n" + body + "}\n"

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name,     params=params,
        body=body, enableCparameters=False)

def add_uniform_Lagrange_interp_3D_to_Cfunction_dict(prefunc=""):
    includes = ["NRPy_basic_defines.h", # // Defines REAL, NinterpGHOSTS, NGHOSTS, & interpstruct. This header is generated in NRPy+.
                "stdio.h", "stdlib.h", "math.h"]
    prefunc += """// The following will result in cache misses just below since the innermost loop is
//     over gf, but it ensures consistency with IDX4() in the rest of NRPy+,
//     and ensures quick access with minimal cache misses when accessed from elsewhere.
#define INTERP_IDX(num_interp_pts,gf,pt) ( (pt) + (num_interp_pts) * (gf) )
"""
    desc = "uniform_Lagrange_interp_3D(): interpolate to an arbitrary list of points list_of_interp_pts_x012"
    name = "uniform_Lagrange_interp_3D"
    params = """const int Nxx_plus_2NGHOSTS0,const int Nxx_plus_2NGHOSTS1,const int Nxx_plus_2NGHOSTS2,  REAL *restrict xx[3],const REAL *dx012_term_inv,

    const REAL *restrict in_gfs,
    const int num_interp_gfs, const int list_of_interp_gfs[num_interp_gfs],
    const int num_interp_pts, const REAL *restrict list_of_interp_pts_x0, const REAL *restrict list_of_interp_pts_x1, const REAL *restrict list_of_interp_pts_x2,

    const int N0,const int N1,const int N2,

    REAL *restrict interp_output"""
    body = r"""REAL *restrict f;
int *restrict x0i_idx;
int *restrict x1j_idx;
int *restrict x2k_idx;
REAL *restrict x0i;
REAL *restrict x1j;
REAL *restrict x2k;
REAL *restrict l0i__times__w0i_inv;
REAL *restrict m1j__times__w1j_inv;
REAL *restrict n2k__times__w2k_inv;
#pragma omp parallel private(f,x0i_idx,x1j_idx,x2k_idx,x0i,x1j,x2k,l0i__times__w0i_inv,m1j__times__w1j_inv,n2k__times__w2k_inv)
{
    f = (REAL *)malloc(sizeof(REAL)*(N0+1)*(N1+1)*(N2+1));
    x0i_idx = (int *)malloc(sizeof(int)*(N0+1));
    x1j_idx = (int *)malloc(sizeof(int)*(N1+1));
    x2k_idx = (int *)malloc(sizeof(int)*(N2+1));
    x0i = (REAL *)malloc(sizeof(REAL)*(N0+1));
    x1j = (REAL *)malloc(sizeof(REAL)*(N1+1));
    x2k = (REAL *)malloc(sizeof(REAL)*(N2+1));
    l0i__times__w0i_inv = (REAL *)malloc(sizeof(REAL)*(N0+1));
    m1j__times__w1j_inv = (REAL *)malloc(sizeof(REAL)*(N1+1));
    n2k__times__w2k_inv = (REAL *)malloc(sizeof(REAL)*(N2+1));

    #pragma omp for
    for(int interp_pt = 0; interp_pt<num_interp_pts; interp_pt++) {
        const REAL curr_xx[3] = { list_of_interp_pts_x0[interp_pt],
                                  list_of_interp_pts_x1[interp_pt],
                                  list_of_interp_pts_x2[interp_pt] };

        // Total cost: 3 * (2 subs + 1 div + 1 add) + 3 * (1 integer sub, div, and add)
        {
            construct_x0i_x1j_x2k_stencil_arrays_3D(
                Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2,xx,
                curr_xx,
                N0,N1,N2, x0i_idx,x1j_idx,x2k_idx);
        }

        // Total cost: (N0+N1+N2) * 1 assignment
        {
            for(int i0=0;i0<=N0;i0++) x0i[i0] = xx[0][x0i_idx[i0]];
            for(int i1=0;i1<=N1;i1++) x1j[i1] = xx[1][x1j_idx[i1]];
            for(int i2=0;i2<=N2;i2++) x2k[i2] = xx[2][x2k_idx[i2]];
        }

        // Total cost: [INTEGER: (N0-1 + N1-1 + N2-1) * (1 assign + 1 mul + 1 sub) + 6 assigns] + 3 typecasts
        //   PLUS (N0 + N1 + N2) * (1 assignment, 1 divide, 1 multiply, and 1 subtract [1 + 3 + 1 + 1 ~ 6 FLOPs])
        Lagrange_interp_coeffs_3D(
                N0,N1,N2, x0i,x1j,x2k, curr_xx,
                l0i__times__w0i_inv,m1j__times__w1j_inv,n2k__times__w2k_inv);

        for(int gf=0;gf<num_interp_gfs;gf++) {
            const int which_gf = list_of_interp_gfs[gf];
            // Total cost: N0N1N2 assignments
            int iii=0;
            for(int i0=0;i0<=N0;i0++) for(int i1=0;i1<=N1;i1++) for(int i2=0;i2<=N2;i2++) {
                  f[iii] = in_gfs[IDX4S(which_gf,x0i_idx[i0],x1j_idx[i1],x2k_idx[i2])];
                  iii++;
                }
            interp_output[INTERP_IDX(num_interp_pts, gf,interp_pt)] =
                (*dx012_term_inv)*
                Lagrange_sum_3D(N0,N1,N2,
                                l0i__times__w0i_inv,
                                m1j__times__w1j_inv,
                                n2k__times__w2k_inv,
                                f);
        }
    }
    free(f);
    free(x0i_idx);
    free(x1j_idx);
    free(x2k_idx);
    free(x0i);
    free(x1j);
    free(x2k);
    free(l0i__times__w0i_inv);
    free(m1j__times__w1j_inv);
    free(n2k__times__w2k_inv);
}
"""
    body = indent_Ccode('    ', body)
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        name=name, params=params,
        body=body, enableCparameters=False)
