from outputC import outputC, add_to_Cfunction_dict


def scalar_field_ID_function_string(Phi, Pi):

    # Set spatial dimension to 3
    # DIM = 3

    # Define lists with lhs/rhs pairs
    rhss = [Phi, Pi]
    lhss = ["in_gfs[IDX4S(PHIGF,i0,i1,i2)]", "in_gfs[IDX4S(PIGF,i0,i1,i2)]"]

    # Get the C function body corresponding to the lhs/rhs lists
    body = outputC(rhss, lhss, filename="returnstring",
                   # outCverbose=False to prevent
                   params="preindent=1,CSE_enable=True,outCverbose=False",
                   # enormous output files.
                   prestring="", poststring="")

    # Add function to dictionary
    desc = "Set up the scalar field initial data at all points on the numerical grid."
    add_to_Cfunction_dict(
        desc=desc,
        name="fields_initial_data",
        params="const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs",
        body=body,
        loopopts="AllPoints,Read_xxs")

    return body
