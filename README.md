# Roc-linalg

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/imclerran/roc-linalg/ci.yaml)

A Roc library for mathematical operations on arbitrary sized vectors and matrices.

```roc
    Matrix,
    # Construction & Conversion
    new,
    from_rows,
    from_columns,
    identity,
    zeros,
    ones,
    concat,
    to_frac,
    # Accessors & Slicing
    dimensions,
    get_rows,
    get_row,
    get_columns,
    get_column,
    element,
    slice,
    minor,
    diagonal,
    diagonal_vector,
    # Row & Column Manipulation
    swap_rows,
    scale_row,
    replace_row,
    drop_row,
    drop_col,
    # Scalar & Elementwise Operations
    scale,
    add_scalar,
    subtract_scalar,
    map_elementwise,
    map2_elementwise,
    multiply_elementwise,
    divide_elementwise,
    # Matrix Arithmetic
    add,
    subtract,
    multiply,
    power,
    transpose,
    multiply_vector,
    # Decompositions, Solvers & Rank
    row_echelon_form,
    row_echelon_form_int,
    reduced_row_echelon_form,
    solve,
    solve_many,
    invert,
    rank,
    # Determinants & Trace
    determinant,
    determinant_int,
    trace,
    # Structural Tests
    is_square,
    is_symmetric,
    is_diagonal,
    is_orthogonal,
    approx_eq,
    # Aggregations & Norms
    sum,
    row_sums,
    column_sums,
    mean,
    row_means,
    column_means,
    frobenius_norm,
```

```roc
    Vector,
    # Constructors & Conversion
    new,
    zeros,
    ones,
    concat,
    to_list,
    to_frac,
    # Accessors
    count,
    element,
    # Scalar Operations
    scale,
    add_scalar,
    subtract_scalar,
    # Vector Arithmetic
    add,
    subtract,
    multiply,
    # Elementwise Operations
    map,
    map2,
    multiply_elementwise,
    divide_elementwise,
    # Aggregations & Averages
    sum,
    mean,
    # Distance, Norms, & Normalization
    distance,
    magnitude,
    l1_norm,
    l2_norm,
    linf_norm,
    normalize,
    # Comparison
    approx_eq,
```