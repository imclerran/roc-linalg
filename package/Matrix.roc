module [
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
]

import rtils.Unsafe exposing [unwrap]

import Vector exposing [Vector]

Matrix a := {
    n_rows : U64,
    n_cols : U64,
    vectors : List (Vector a),
}
    implements [Eq, Inspect]

#==============================================================================
# Construction & Conversion
#==============================================================================

## Create a new matrix from a list of lists (rows) of numbers.
new : List (List (Num a)) -> Result (Matrix a) [RowLengthMismatch, MatrixCannotBeEmpty]
new = |ll|
    vectors = List.map(ll, Vector.new)
    n_rows = List.len(ll)
    vector_lengths = List.map(vectors, Vector.count)
    if List.min(vector_lengths) != List.max(vector_lengths) then
        Err(RowLengthMismatch)
    else
        n_cols = List.first(vector_lengths) ? |ListWasEmpty| MatrixCannotBeEmpty
        @Matrix(
            {
                n_rows,
                n_cols,
                vectors,
            },
        )
        |> Ok

## Create a matrix from a list of row vectors.
from_rows : List (Vector a) -> Result (Matrix a) [RowLengthMismatch, MatrixCannotBeEmpty]
from_rows = |vectors|
    min_len = List.map(vectors, Vector.count) |> List.min |> Result.with_default(0)
    max_len = List.map(vectors, Vector.count) |> List.max |> Result.with_default(0)
    if min_len != max_len then
        Err(RowLengthMismatch)
    else if min_len < 1 then
        Err(MatrixCannotBeEmpty)
    else
        @Matrix(
            {
                n_rows: List.len(vectors),
                n_cols: min_len,
                vectors,
            },
        )
        |> Ok

## Create a matrix from a list of column vectors.
from_columns : List (Vector a) -> Result (Matrix a) [ColumnLengthMismatch, MatrixCannotBeEmpty]
from_columns = |cols|
    min_len = List.map(cols, Vector.count) |> List.min |> Result.with_default(0)
    max_len = List.map(cols, Vector.count) |> List.max |> Result.with_default(0)
    if min_len != max_len then
        Err(ColumnLengthMismatch)
    else if min_len < 1 then
        Err(MatrixCannotBeEmpty)
    else
        vectors =
            List.range({ start: At(0), end: Length(min_len) })
            |> List.walk(
                [],
                |vs, r|
                    List.append(
                        vs,
                        List.map(
                            cols,
                            |col|
                                Vector.element(col, r) |> unwrap("row index is always inbounds"),
                        )
                        |> Vector.new,
                    ),
            )
        @Matrix(
            {
                n_cols: List.len(cols),
                n_rows: min_len,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) |> unwrap("failed to construct matrix mx1")
    cs = get_columns(mx1)
    mx2 = from_columns(cs)
    Ok(mx1) == mx2

## Create an identity matrix of the given size (n x n).
identity : U64 -> Matrix a
identity = |n|
    vectors =
        List.range({ start: At(0), end: Length(n) })
        |> List.map(
            |i|
                List.range({ start: At(0), end: Length(n) })
                |> List.map(|j| if i == j then 1 else 0),
        )
        |> List.map(Vector.new)
    @Matrix(
        {
            n_rows: n,
            n_cols: n,
            vectors,
        },
    )

expect new([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == Ok(identity(3))

## Create a matrix of the given dimensions filled with zeros.
zeros : U64, U64 -> Matrix a
zeros = |n_rows, n_cols|
    vectors = List.repeat(Vector.zeros(n_cols), n_rows)
    @Matrix(
        {
            n_rows,
            n_cols,
            vectors,
        },
    )

## Create a matrix of the given dimensions filled with ones.
ones : U64, U64 -> Matrix a
ones = |n_rows, n_cols|
    vectors = List.repeat(Vector.ones(n_cols), n_rows)
    @Matrix(
        {
            n_rows,
            n_cols,
            vectors,
        },
    )

## Concatenate two matrices horizontally.
concat : Matrix a, Matrix a -> Result (Matrix a) [DimensionMismatch]
concat = |@Matrix(mx1), @Matrix(mx2)|
    if mx1.n_rows != mx2.n_rows then
        Err(DimensionMismatch)
    else
        vectors = List.map2(mx1.vectors, mx2.vectors, Vector.concat)
        @Matrix(
            {
                n_rows: mx1.n_rows,
                n_cols: mx1.n_cols + mx2.n_cols,
                vectors,
            },
        )
        |> Ok

## Convert a matrix of numbers to a matrix of floating-point numbers.
to_frac : Matrix * -> Matrix (FloatingPoint *)
to_frac = |@Matrix(mx)|
    vectors = List.map(mx.vectors, Vector.to_frac)
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

#==============================================================================
# Accessors & Slicing
#==============================================================================

## Get the dimensions of the matrix as a record with row and column counts.
dimensions : Matrix a -> { rows : U64, cols : U64 }
dimensions = |@Matrix(mx)| { rows: mx.n_rows, cols: mx.n_cols }

## Get the rows of the matrix as a list of vectors.
get_rows : Matrix a -> List (Vector a)
get_rows = |@Matrix(mx)| mx.vectors

## Get the columns of the matrix as a list of vectors.
get_columns : Matrix a -> List (Vector a)
get_columns = |@Matrix(mx)|
    res =
        List.range({ start: At 0, end: Length(mx.n_cols) })
        |> List.map_try(
            |c|
                get_column(@Matrix(mx), c),
        )
    when res is
        Ok(cs) -> cs
        _ -> crash "index is always in range"

## Get the row at the specified index.
get_row : Matrix a, U64 -> Result (Vector a) [OutOfBounds]
get_row = |@Matrix(mx), index| List.get(mx.vectors, index)

## Get the column at the specified index.
get_column : Matrix a, U64 -> Result (Vector a) [OutOfBounds]
get_column = |@Matrix(mx), index| List.map_try(mx.vectors, |r| Vector.element(r, index)) |> Result.map_ok(Vector.new)

## Get the element at the specified row and column indices.
element : Matrix a, U64, U64 -> Result (Num a) [OutOfBounds]
element = |@Matrix(mx), r, c| List.get(mx.vectors, r)? |> Vector.element(c)

Range : {
    start : [At U64, After U64],
    end : [At U64, Before U64, Length U64],
}

## Get a slice of the matrix.
## ```
## Range : {
##     start : [At U64, After U64],
##     end : [At U64, Before U64, Length U64],
## }
## ```
slice : Matrix a, Range, Range -> Result (Matrix a) [InvalidRowRange, InvalidColumnRange, MatrixCannotBeEmpty]
slice = |@Matrix(mx), rows, cols|
    to_start_and_len = |range|
        start =
            when range.start is
                At(n) -> n
                After(n) -> n + 1
        when range.end is
            At(n) ->
                Ok({ start, len: n - start })

            Before(n) if n > (start + 1) ->
                Ok({ start, len: (n - start) - 1 })

            Before(_) ->
                Err(LengthLessThanOne)

            Length(n) ->
                Ok({ start, len: n })

    range_r = to_start_and_len(rows) ? |LengthLessThanOne| InvalidRowRange
    range_c = to_start_and_len(cols) ? |LengthLessThanOne| InvalidColumnRange

    vectors =
        List.sublist(mx.vectors, range_r)
        |> List.map(
            |r|
                Vector.to_list(r)
                |> List.sublist(range_c)
                |> Vector.new,
        )

    @Matrix(
        {
            n_rows: mx.n_rows - range_r.len,
            n_cols: mx.n_cols - range_c.len,
            vectors,
        },
    )
    |> Ok

## Get the minor matrix by removing the specified row and column.
minor : Matrix a, U64, U64 -> Result (Matrix a) [MatrixCannotBeEmpty]
minor = |mx, r, c| mx |> drop_row(r)? |> drop_col(c)

expect
    mx1 = new([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) |> unwrap("failed to contruct matrix mx1")
    mx2 = minor(mx1, 1, 1)
    mx2 == new([[1, 3], [7, 9]])

## Get the diagonal matrix formed from the diagonal elements of the input matrix.
diagonal : Matrix a -> Result (Matrix a) [MatrixMustBeSquare]
diagonal = |@Matrix(mx)|
    if !is_square(@Matrix(mx)) then
        Err(MatrixMustBeSquare)
    else
        vectors =
            List.walk_with_index(
                mx.vectors,
                [],
                |vecs, vec, i|
                    v =
                        Vector.to_list(vec)
                        |> List.walk_with_index(
                            [],
                            |elems, el, j|
                                if i == j then
                                    List.append(elems, el)
                                else
                                    List.append(elems, 0),
                        )
                        |> Vector.new
                    List.append(vecs, v),
            )
        @Matrix(
            {
                n_rows: mx.n_rows,
                n_cols: mx.n_cols,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[1, 1], [1, 1]]) |> unwrap("failed to construct matrix mx1")
    mx2 = diagonal(mx1)
    mx3 = new([[1, 0], [0, 1]])
    mx2 == mx3

## Get the diagonal elements of the matrix as a vector.
diagonal_vector : Matrix a -> Result (Vector a) [MatrixMustBeSquare]
diagonal_vector = |@Matrix(mx)|
    if !is_square(@Matrix(mx)) then
        Err(MatrixMustBeSquare)
    else
        vals =
            List.range({ start: At(0), end: Length(mx.n_rows) })
            |> List.map(
                |i|
                    @Matrix(mx)
                    |> element(i, i)
                    |> unwrap("diagonal element always inbounds"),
            )
        Vector.new(vals) |> Ok

#==============================================================================
# Row & Column Manipulation
#==============================================================================

## Swap two rows in the matrix.
swap_rows : Matrix a, U64, U64 -> Matrix a
swap_rows = |@Matrix(mx), r1, r2|
    vectors = List.swap(mx.vectors, r1, r2)
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

## Multiply a row by a scalar.
scale_row : Matrix a, U64, Num a -> Matrix a
scale_row = |@Matrix(mx), r, scalar|
    vectors =
        List.map_with_index(
            mx.vectors,
            |vec, i|
                if i == r then
                    Vector.scale(vec, scalar)
                else
                    vec,
        )
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

## Replace a row with a new vector.
replace_row : Matrix a, U64, Vector a -> Matrix a
replace_row = |@Matrix(mx), r, new_vec|
    vectors =
        List.map_with_index(
            mx.vectors,
            |vec, i|
                if i == r then
                    new_vec
                else
                    vec,
        )
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

## Drop a row from the matrix.
drop_row : Matrix a, U64 -> Result (Matrix a) [MatrixCannotBeEmpty]
drop_row = |@Matrix(mx), r|
    vectors = List.drop_at(mx.vectors, r)
    n_rows = List.len(vectors)
    if n_rows == 0 then
        Err(MatrixCannotBeEmpty)
    else
        @Matrix(
            {
                n_rows: List.len(vectors),
                n_cols: mx.n_cols,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[1, 1], [2, 2], [3, 3]]) |> unwrap("failed to construct matrix mx1")
    mx2 = drop_row(mx1, 1)
    mx2 == new([[1, 1], [3, 3]])

## Drop a column from the matrix.
drop_col : Matrix a, U64 -> Result (Matrix a) [MatrixCannotBeEmpty]
drop_col = |@Matrix(mx), c|
    vectors =
        List.map(
            mx.vectors,
            |vec|
                Vector.to_list(vec)
                |> List.walk_with_index(
                    [],
                    |vals, v, i|
                        if i == c then
                            vals
                        else
                            List.append(vals, v),
                )
                |> Vector.new,
        )
    n_cols =
        List.first(vectors)
        |> Result.map_ok(Vector.to_list)
        |> Result.map_ok(List.len)
        |> Result.map_err(|ListWasEmpty| MatrixCannotBeEmpty)?

    if n_cols == 0 then
        Err(MatrixCannotBeEmpty)
    else
        @Matrix(
            {
                n_cols,
                n_rows: mx.n_rows,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[1, 2, 3], [1, 2, 3]]) |> unwrap("failed to construct matrix mx1")
    mx2 = drop_col(mx1, 1)
    mx2 == new([[1, 3], [1, 3]])

#==============================================================================
# Scalar & Elementwise Operations
#==============================================================================

## Multiply the matrix by a scalar.
scale : Matrix a, Num a -> Matrix a
scale = |@Matrix(mx), n|
    vectors = List.map(mx.vectors, |v| Vector.scale(v, n))
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

## Add a scalar to each element of the matrix.
add_scalar : Matrix a, Num a -> Matrix a
add_scalar = |@Matrix(mx), n|
    vectors = List.map(mx.vectors, |v| Vector.add_scalar(v, n))
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

## Subtract a scalar from each element of the matrix.
subtract_scalar : Matrix a, Num a -> Matrix a
subtract_scalar = |@Matrix(mx), n|
    vectors = List.map(mx.vectors, |v| Vector.subtract_scalar(v, n))
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

## Multiply two matrices elementwise.
multiply_elementwise : Matrix a, Matrix a -> Result (Matrix a) [DimensionMismatch]
multiply_elementwise = |mx1, mx2| map2_elementwise(mx1, mx2, Num.mul)

## Divide two matrices elementwise.
divide_elementwise : Matrix (FloatingPoint a), Matrix (FloatingPoint a) -> Result (Matrix (FloatingPoint a)) [DimensionMismatch]
divide_elementwise = |mx1, mx2| map2_elementwise(mx1, mx2, Num.div)

## Element-wise map of a matrix with a unary function.
map_elementwise : Matrix a, (Num a -> Num a) -> Matrix a
map_elementwise = |@Matrix(mx), transform|
    vectors =
        List.map(
            mx.vectors,
            |v|
                Vector.to_list(v)
                |> List.map(transform)
                |> Vector.new,
        )
    @Matrix(
        { mx &
            vectors,
        },
    )

## Element-wise map of two matrices with a binary function.
map2_elementwise : Matrix a, Matrix b, (Num a, Num b -> Num c) -> Result (Matrix c) [DimensionMismatch]
map2_elementwise = |@Matrix(mx1), @Matrix(mx2), combine|
    if dimensions(@Matrix(mx1)) != dimensions(@Matrix(mx2)) then
        Err(DimensionMismatch)
    else
        vectors =
            List.map2(
                mx1.vectors,
                mx2.vectors,
                |v1, v2|
                    List.map2(
                        Vector.to_list(v1),
                        Vector.to_list(v2),
                        combine,
                    )
                    |> Vector.new,
            )
        @Matrix(
            {
                n_rows: mx1.n_rows,
                n_cols: mx1.n_cols,
                vectors,
            },
        )
        |> Ok

#==============================================================================
# Matrix Arithmetic
#==============================================================================

## Add two matrices.
add : Matrix a, Matrix a -> Result (Matrix a) [DimensionMismatch]
add = |@Matrix(mx1), @Matrix(mx2)|
    if mx1.n_rows != mx2.n_rows or mx1.n_cols != mx2.n_cols then
        Err(DimensionMismatch)
    else
        vectors =
            List.map2(mx1.vectors, mx2.vectors, |r1, r2| Vector.add(r1, r2))
            |> List.map_try(|v| v)
            |> Result.map_err(|VectorLengthMismatch| DimensionMismatch)?
        @Matrix(
            {
                n_rows: mx1.n_rows,
                n_cols: mx1.n_cols,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[1, 2], [3, 4]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[5, 6], [7, 8]]) |> unwrap("failed to construct matrix mx2")
    mx3 = new([[6, 8], [10, 12]])
    mx3 == add(mx1, mx2)

## Subtract two matrices.
subtract : Matrix a, Matrix a -> Result (Matrix a) [DimensionMismatch]
subtract = |@Matrix(mx1), @Matrix(mx2)|
    if mx1.n_rows != mx2.n_rows or mx1.n_cols != mx2.n_cols then
        Err(DimensionMismatch)
    else
        vectors =
            List.map2(mx1.vectors, mx2.vectors, |r1, r2| Vector.subtract(r1, r2))
            |> List.map_try(|v| v)
            |> Result.map_err(|VectorLengthMismatch| DimensionMismatch)?
        @Matrix(
            {
                n_rows: mx1.n_rows,
                n_cols: mx1.n_cols,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[5, 6], [7, 8]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[1, 2], [3, 4]]) |> unwrap("failed to construct matrix mx2")
    mx3 = new([[4, 4], [4, 4]])
    mx3 == subtract(mx1, mx2)

## Perform matrix multiplication on two matrices.
multiply : Matrix a, Matrix a -> Result (Matrix a) [DimensionMismatch]
multiply = |@Matrix(mx1), @Matrix(mx2)|
    if mx1.n_cols != mx2.n_rows then
        Err(DimensionMismatch)
    else
        vectors =
            List.map_try(mx1.vectors, |r| List.map_try(get_columns(@Matrix(mx2)), |c| Vector.multiply(r, c)))
            |> Result.map_err(|VectorLengthMismatch| DimensionMismatch)?
            |> List.map(Vector.new)
        @Matrix(
            {
                n_rows: mx1.n_rows,
                n_cols: mx2.n_cols,
                vectors,
            },
        )
        |> Ok

expect
    mx1 = new([[0, 1, 2], [-3, 2, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[-1, 2], [2, 1], [3, 0]]) |> unwrap("failed to construct matrix mx2")
    mx3 = new([[8, 1], [7, -4]])
    mx3 == multiply(mx1, mx2)

## Raise a matrix to an integer power.
power : Matrix a, U64 -> Result (Matrix a) [ZeroPowerOnNonSquareMatrix, DimensionMismatch]
power = |@Matrix(mx), p|
    when p is
        0 ->
            if mx.n_rows != mx.n_cols then
                Err(ZeroPowerOnNonSquareMatrix)
            else
                identity(mx.n_rows) |> Ok

        1 -> @Matrix(mx) |> Ok
        _ -> power(multiply(@Matrix(mx), @Matrix(mx))?, p - 1)

expect
    mx1 = new([[1, 2], [1, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = identity(2)
    Ok(mx2) == power(mx1, 0)

expect
    mx1 = new([[1, 2], [1, 0]]) |> unwrap("failed to construct matrix mx1")
    Ok(mx1) == power(mx1, 1)

expect
    mx1 = new([[1, 2], [1, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[3, 2], [1, 2]])
    mx2 == power(mx1, 2)

expect
    mx0 = identity(2)
    mx1 = new([[1, 2], [1, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = power(mx1, 2)
    mx2 == multiply(mx1, mx1) |> Result.try(|mx| multiply(mx, mx0))

## Transpose the matrix.
transpose : Matrix a -> Matrix a
transpose = |@Matrix(mx)|
    @Matrix(
        {
            n_rows: mx.n_cols,
            n_cols: mx.n_rows,
            vectors: get_columns(@Matrix(mx)),
        },
    )

expect
    mx1 = new([[1, 2, 3], [4, 5, 6]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[1, 4], [2, 5], [3, 6]]) |> unwrap("failed to construct matrix mx2")
    mx2 == transpose(mx1)

## Multiply a matrix by a vector.
multiply_vector : Matrix a, Vector a -> Result (Vector a) [DimensionMismatch]
multiply_vector = |@Matrix(mx), vec|
    List.map_try(mx.vectors, |v| Vector.multiply(vec, v))
    |> Result.map_err(|_| DimensionMismatch)
    |> Result.map_ok(Vector.new)

#==============================================================================
# Decompositions, Solvers & Rank
#==============================================================================

## Compute the row echelon form of a matrix, returning a matrix of floating-point numbers.  
row_echelon_form : Matrix a -> Matrix (FloatingPoint b)
row_echelon_form = |@Matrix(mx)|
    get_pivot_idx = |v| Vector.to_list(v) |> List.find_first_index(|x| !Num.is_approx_eq(x, 0, { atol: 0.000_000_1 }))

    compare_rows = |v1, v2|
        pivot1 = get_pivot_idx(v1) |> Result.with_default(mx.n_cols)
        pivot2 = get_pivot_idx(v2) |> Result.with_default(mx.n_cols)
        if pivot1 < pivot2 then
            LT
        else if pivot1 > pivot2 then
            GT
        else
            EQ

    sorted_vectors =
        List.map(mx.vectors, Vector.to_frac)
        |> List.sort_with(compare_rows)

    vectors =
        List.walk_with_index(
            sorted_vectors,
            sorted_vectors,
            |acc, _, r|
                pivot_vec = List.get(acc, r) |> unwrap("row vector always inbounds")
                idx_res = get_pivot_idx(pivot_vec)
                when idx_res is
                    Ok(pivot_idx) ->
                        pivot_val =
                            Vector.element(pivot_vec, pivot_idx)
                            |> unwrap("pivot index is always inbounds")
                        List.walk_with_index(
                            acc,
                            [],
                            |acc2, vec2, r2|
                                if r2 > r then
                                    val2 =
                                        Vector.element(vec2, pivot_idx)
                                        |> unwrap("pivot index is always inbounds")
                                    scalar = val2 / pivot_val
                                    transformed =
                                        Vector.subtract(vec2, Vector.scale(pivot_vec, scalar))
                                        |> unwrap("vectors are always the same length")
                                    List.append(acc2, transformed)
                                else
                                    List.append(acc2, vec2),
                        )

                    Err(_no_pivot_in_row) -> acc,
        )
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

expect
    mx1 = new([[3, 1, 2], [2, 1, 3], [0, 0, 0], [0, 1, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[3, 1, 2], [0, 1 / 3, 5 / 3], [0, 0, -5], [0, 0, 0]]) |> unwrap("failed to construct matrix mx2")
    mx3 = row_echelon_form(mx1)
    mx2 |> approx_eq(mx3)

## Compute the row echelon form of an integer matrix, returning a matrix of integers.
row_echelon_form_int : Matrix (Integer a) -> Matrix (Integer a)
row_echelon_form_int = |@Matrix(mx)|
    least_common_multiple = |n1, n2| n1 * n2 // greatest_common_divisor(n1, n2)

    greatest_common_divisor = |n1, n2| if Num.is_zero(n1) then n2 else greatest_common_divisor(n2 % n1, n1)

    get_pivot_idx = |v| Vector.to_list(v) |> List.find_first_index(|x| x != 0)

    compare_rows = |v1, v2|
        pivot1 = get_pivot_idx(v1) |> Result.with_default(mx.n_cols)
        pivot2 = get_pivot_idx(v2) |> Result.with_default(mx.n_cols)
        if pivot1 < pivot2 then
            LT
        else if pivot1 > pivot2 then
            GT
        else
            EQ

    sorted_vectors = List.sort_with(mx.vectors, compare_rows)

    vectors =
        List.walk_with_index(
            sorted_vectors,
            sorted_vectors,
            |acc, _, r|
                pivot_vec = List.get(acc, r) |> unwrap("row vector always inbounds")
                idx_res = get_pivot_idx(pivot_vec)
                when idx_res is
                    Ok(pivot_idx) ->
                        pivot_val =
                            Vector.element(pivot_vec, pivot_idx)
                            |> unwrap("pivot index is always inbounds")
                        List.walk_with_index(
                            acc,
                            [],
                            |acc2, vec2, r2|
                                if r2 > r then
                                    val2 =
                                        Vector.element(vec2, pivot_idx)
                                        |> unwrap("pivot index is always inbounds")
                                    if val2 == 0 then
                                        List.append(acc2, vec2)
                                    else
                                        lcm = least_common_multiple(Num.abs(pivot_val), Num.abs(val2))
                                        scalar1 = lcm // Num.abs(val2)
                                        scalar2 = lcm // Num.abs(pivot_val)
                                        if (pivot_val > 0) == (val2 > 0) then
                                            transformed =
                                                Vector.subtract(
                                                    Vector.scale(vec2, scalar1),
                                                    Vector.scale(pivot_vec, scalar2),
                                                )
                                                |> unwrap("vectors are always the same length")
                                            List.append(acc2, transformed)
                                        else
                                            transformed =
                                                Vector.add(
                                                    Vector.scale(vec2, scalar1),
                                                    Vector.scale(pivot_vec, scalar2),
                                                )
                                                |> unwrap("vectors are always the same length")
                                            List.append(acc2, transformed)
                                else
                                    List.append(acc2, vec2),
                        )

                    Err(_no_pivot_in_row) -> acc,
        )
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

expect
    mx1 = new([[3, 1, 2], [2, 1, 3], [0, 0, 0], [0, 1, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[3, 1, 2], [0, 1, 5], [0, 0, -5], [0, 0, 0]]) |> unwrap("failed to construct matrix mx2")
    mx3 = row_echelon_form_int(mx1)
    mx2 == mx3

## Compute the reduced row echelon form of a matrix, returning a matrix of floating-point numbers.
reduced_row_echelon_form : Matrix a -> Matrix (FloatingPoint b)
reduced_row_echelon_form = |@Matrix(mx)|
    get_pivot_idx = |v| Vector.to_list(v) |> List.find_first_index(|x| !Num.is_approx_eq(x, 0, { atol: 0.000_000_1 }))

    compare_rows = |v1, v2|
        pivot1 = get_pivot_idx(v1) |> Result.with_default(mx.n_cols)
        pivot2 = get_pivot_idx(v2) |> Result.with_default(mx.n_cols)
        if pivot1 < pivot2 then
            LT
        else if pivot1 > pivot2 then
            GT
        else
            EQ

    sorted_vectors =
        List.map(mx.vectors, Vector.to_frac)
        |> List.sort_with(compare_rows)

    vectors =
        List.walk_with_index(
            sorted_vectors,
            sorted_vectors,
            |acc, _, r|
                pivot_vec = List.get(acc, r) |> unwrap("row vector always inbounds")
                idx_res = get_pivot_idx(pivot_vec)
                when idx_res is
                    Ok(pivot_idx) ->
                        pivot_val =
                            Vector.element(pivot_vec, pivot_idx)
                            |> unwrap("pivot index is always inbounds")
                        normalized_pivot_row = Vector.scale(pivot_vec, 1 / pivot_val)
                        List.walk_with_index(
                            acc,
                            [],
                            |acc2, vec2, r2|
                                if r2 != r then
                                    val2 =
                                        Vector.element(vec2, pivot_idx)
                                        |> unwrap("pivot index is always inbounds")
                                    scalar = val2
                                    transformed =
                                        Vector.subtract(vec2, Vector.scale(normalized_pivot_row, scalar))
                                        |> unwrap("vectors are always the same length")
                                    List.append(acc2, transformed)
                                else
                                    List.append(acc2, normalized_pivot_row),
                        )

                    Err(_no_pivot_in_row) -> acc,
        )
    @Matrix(
        {
            n_rows: mx.n_rows,
            n_cols: mx.n_cols,
            vectors,
        },
    )

expect
    mx1 = new([[3,1,2], [2,1,3]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[1,0, -1], [0,1, 5]]) |> unwrap("failed to construct matrix mx2")
    mx3 = reduced_row_echelon_form(mx1)
    mx2 |> approx_eq(mx3)

## Solve the linear system Ax = b for x, where A is the matrix and b is the vector.
solve : Matrix a, Vector a -> Result (Vector (FloatingPoint b)) [MatrixMustBeSquare, DimensionMismatch, DeterminantIsZero]
solve = |mx, b|
    inv = invert(mx)?
    multiply_vector(inv, Vector.to_frac(b))
    |> Result.map_ok(|vec| Vector.to_frac(vec))

## Solve the linear system AX = B for X, where A, B and X are matrices.
solve_many : Matrix a, Matrix a -> Result (Matrix (FloatingPoint b)) [MatrixMustBeSquare, DimensionMismatch, DeterminantIsZero]
solve_many = |mx, b|
    inv = invert(mx)?
    multiply(inv, to_frac(b))

## Invert a square matrix.
invert : Matrix * -> Result (Matrix (FloatingPoint *)) [MatrixMustBeSquare, DeterminantIsZero]
invert = |mx|
    det = determinant(mx)?
    if det == 0 then
        Err(DeterminantIsZero)
    else
        dim = dimensions(mx)
        left = to_frac(mx)
        right = identity(dim.rows) |> to_frac
        aug = concat(left, right) |> unwrap("identity matrix has same height as initial matrix")

        # Perform Gauss-Jordan elimination to reduce [A|I] -> [I|A^{-1}]
        final_rows =
            List.walk_with_index(
                get_rows(aug),
                get_rows(aug),
                |rows_acc, _, c|
                    # Find a pivot row r >= c with non-zero at column c
                    r =
                        List.range({ start: At(c), end: Length(dim.rows - c) })
                        |> List.walk(
                            c,
                            |acc, i|
                                if acc != c then
                                    acc
                                else
                                    vec_i = List.get(rows_acc, i) |> unwrap("row index is always inbounds")
                                    x = Vector.element(vec_i, c) |> unwrap("column index is always inbounds")
                                    if Num.is_approx_eq(x, Num.to_frac(0), {}) then
                                        acc
                                    else
                                        i,
                        )

                    rows_swapped = if r == c then rows_acc else List.swap(rows_acc, c, r)

                    # Normalize pivot row so pivot becomes 1
                    pivot_val =
                        List.get(rows_swapped, c)
                        |> unwrap("pivot row index is always inbounds")
                        |> Vector.element(c)
                        |> unwrap("pivot column index is always inbounds")
                    pivot_row =
                        List.get(rows_swapped, c)
                        |> unwrap("pivot row index is always inbounds")
                        |> Vector.scale(1 / pivot_val)

                    # Eliminate column c from all other rows
                    List.map_with_index(
                        rows_swapped,
                        |vec, i|
                            if i == c then
                                pivot_row
                            else
                                scalar =
                                    Vector.element(vec, c)
                                    |> unwrap("element index always inbounds")
                                Vector.subtract(vec, Vector.scale(pivot_row, scalar))
                                |> unwrap("vectors are always the same length"),
                    ),
            )

        inv_vectors =
            List.map(
                final_rows,
                |vec|
                    vals = Vector.to_list(vec)
                    right_vals = List.drop_first(vals, dim.cols)
                    Vector.new(right_vals),
            )

        @Matrix(
            {
                n_rows: dim.rows,
                n_cols: dim.cols,
                vectors: inv_vectors,
            },
        )
        |> Ok

## Compute the rank of the matrix.
rank : Matrix a -> U64
rank = |mx|
    has_pivot = |v|
        Vector.to_list(v)
        |> List.find_first_index(|x| !Num.is_approx_eq(x, 0, { atol: 0.000_000_1 }))
        |> Result.is_ok

    reduced = to_frac(mx) |> row_echelon_form

    List.walk_until(
        get_rows(reduced),
        0,
        |count, vec| if has_pivot(vec) then Continue(count + 1) else Break(count),
    )

#==============================================================================
# Determinants & Trace
#==============================================================================

## Compute the determinant of a square matrix.
determinant : Matrix * -> Result (Frac *) [MatrixMustBeSquare]
determinant = |@Matrix(mx)|
    if mx.n_cols != mx.n_rows then
        Err(MatrixMustBeSquare)
    else if mx.n_cols == 1 then
        @Matrix(mx) |> element(0, 0) |> unwrap("idx 0,0 inbounds for 1x1 matrix") |> Num.to_frac |> Ok
    else if mx.n_cols == 2 then
        a = @Matrix(mx) |> element(0, 0) |> unwrap("idx 0,0 inbounds for 2x2 matrix") |> Num.to_frac
        b = @Matrix(mx) |> element(0, 1) |> unwrap("idx 0,1 inbounds for 2x2 matrix") |> Num.to_frac
        c = @Matrix(mx) |> element(1, 0) |> unwrap("idx 1,0 inbounds for 2x2 matrix") |> Num.to_frac
        d = @Matrix(mx) |> element(1, 1) |> unwrap("idx 1,1 inbounds for 2x2 matrix") |> Num.to_frac
        (a * d) - (b * c) |> Ok
    else
        List.range({ start: At(0), end: Length(mx.n_cols) })
        |> List.walk(
            0,
            |det, i|
                m = minor(@Matrix(mx), 0, i) |> unwrap("minor from matrix of size > 2x2 is not empty")
                sign = Num.pow(-1, Num.to_frac(i))
                elem = @Matrix(mx) |> element(0, i) |> unwrap("index from 0 to length(len) always inbounds") |> Num.to_frac
                det + sign * elem * (determinant(m) |> unwrap("matrix is square")),
        )
        |> Ok

## Compute the determinant of a square integer matrix.
determinant_int : Matrix (Integer a) -> Result I128 [MatrixMustBeSquare]
determinant_int = |@Matrix(mx)|
    if mx.n_cols != mx.n_rows then
        Err(MatrixMustBeSquare)
    else if mx.n_cols == 1 then
        @Matrix(mx) |> element(0, 0) |> unwrap("idx 0,0 inbounds for 1x1 matrix") |> Num.to_i128 |> Ok
    else if mx.n_cols == 2 then
        a = @Matrix(mx) |> element(0, 0) |> unwrap("idx 0,0 inbounds for 2x2 matrix") |> Num.to_i128
        b = @Matrix(mx) |> element(0, 1) |> unwrap("idx 0,1 inbounds for 2x2 matrix") |> Num.to_i128
        c = @Matrix(mx) |> element(1, 0) |> unwrap("idx 1,0 inbounds for 2x2 matrix") |> Num.to_i128
        d = @Matrix(mx) |> element(1, 1) |> unwrap("idx 1,1 inbounds for 2x2 matrix") |> Num.to_i128
        (a * d) - (b * c) |> Ok
    else
        List.range({ start: At(0), end: Length(mx.n_cols) })
        |> List.walk(
            0,
            |det, i|
                m = minor(@Matrix(mx), 0, i) |> unwrap("minor from matrix of size > 2x2 is not empty")
                sign = Num.pow_int(-1, Num.int_cast(i)) |> Num.to_i128
                elem = @Matrix(mx) |> element(0, i) |> unwrap("index from 0 to length(len) always inbounds") |> Num.to_i128
                det + sign * elem * (determinant_int(m) |> unwrap("matrix is square")),
        )
        |> Ok

expect
    mx1 = new([[1, 2, 3], [3, 2, 1], [2, 1, 3]]) |> unwrap("failed to construct matrix mx1")
    d = determinant_int(mx1)
    Ok(-12) == d

## Compute the trace of a square matrix.
trace : Matrix a -> Result (Num a) [MatrixMustBeSquare]
trace = |mx| diagonal_vector(mx) |> Result.map_ok(|vec| Vector.sum(vec))

##==============================================================================
# Structural Tests
#==============================================================================

## Check if the matrix is square.
is_square : Matrix a -> Bool
is_square = |@Matrix(mx)| mx.n_rows == mx.n_cols

## Check if the matrix is symmetric.
is_symmetric : Matrix a -> Bool
is_symmetric = |mx| is_square(mx) and mx == transpose(mx)

expect
    mx1 = new([[1, 2, 3], [2, 4, 5], [3, 5, 6]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[1, 2, 3], [3, 4, 5], [3, 5, 6]]) |> unwrap("failed to construct matrix mx2")
    is_symmetric(mx1) and !is_symmetric(mx2)

## Check if the matrix is diagonal.
is_diagonal : Matrix a -> Bool
is_diagonal = |@Matrix(mx)|
    if !is_square(@Matrix(mx)) then
        Bool.false
    else
        List.walk_with_index_until(
            mx.vectors,
            Bool.true,
            |_, vec, i|
                ok =
                    List.walk_with_index_until(
                        Vector.to_list(vec),
                        Bool.true,
                        |_, v, j|
                            if i == j or v == 0 then
                                Continue(Bool.true)
                            else
                                Break(Bool.false),
                    )
                if ok then
                    Continue(ok)
                else
                    Break(ok),
        )

expect
    mx1 = new([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[1, 0, 1], [0, 1, 1], [0, 0, 1]]) |> unwrap("failed to construct matrix mx2")
    is_diagonal(mx1) and !is_diagonal(mx2)

## Check if the matrix is orthogonal.
is_orthogonal : Matrix a -> Bool
is_orthogonal = |mx|
    if !is_square(mx) then
        Bool.false
    else
        n = dimensions(mx).rows
        when multiply(transpose(mx), mx) is
            Ok(prod) -> prod == identity(n)
            Err(_) -> Bool.false

expect
    mx1 = new([[0, 1], [-1, 0]]) |> unwrap("failed to construct matrix mx1")
    mx2 = new([[0, 2], [-2, 0]]) |> unwrap("failed to construct matrix mx2")
    is_orthogonal(mx1) and !is_orthogonal(mx2)

## Check if two matrices are approximately equal (for floating-point matrices).
approx_eq : Matrix (FloatingPoint a), Matrix (FloatingPoint a) -> Bool
approx_eq = |@Matrix(mx1), @Matrix(mx2)|
    List.map2(
        mx1.vectors,
        mx2.vectors,
        |v1, v2|
            Vector.approx_eq(v1, v2),
    )
    |> List.all(|b| b)

expect
    mx1 = new([[1, 2, 3], [3, 2, 1], [2, 1, 3]]) |> unwrap("failed to construct matrix mx1")
    mx2 = invert(mx1) |> unwrap("invert produced error message")
    mx3 = new([[-5 / 12, 1 / 4, 1 / 3], [7 / 12, 1 / 4, -2 / 3], [1 / 12, -1 / 4, 1 / 3]]) |> unwrap("failed to construct matrix mx3")
    mx2 |> approx_eq(mx3)

#==============================================================================
# Aggregations & Norms
#==============================================================================

## Compute the sum of all elements in the matrix.
sum : Matrix a -> Num a
sum = |@Matrix(mx)|
    List.map(mx.vectors, Vector.sum)
    |> List.sum

## Compute the sum of each row as a vector.
row_sums : Matrix a -> Vector a
row_sums = |@Matrix(mx)| List.map(mx.vectors, Vector.sum) |> Vector.new

## Compute the sum of each column as a vector.
column_sums : Matrix a -> Vector a
column_sums = |mx| get_columns(mx) |> List.map(Vector.sum) |> Vector.new

## Compute the mean of all elements in the matrix.
mean : Matrix a -> Frac b
mean = |mx|
    total_sum = sum(mx) |> Num.to_frac
    dim = dimensions(mx)
    total_elements = Num.to_frac(dim.rows * dim.cols)
    Num.div(total_sum, total_elements)

## Compute the mean of each row as a vector.
row_means : Matrix a -> Vector (FloatingPoint b)
row_means = |@Matrix(mx)| List.map(mx.vectors, Vector.mean) |> Vector.new

## Compute the mean of each column as a vector.
column_means : Matrix a -> Vector (FloatingPoint b)
column_means = |mx| get_columns(mx) |> List.map(Vector.mean) |> Vector.new

## Compute the Frobenius norm of the matrix.
frobenius_norm : Matrix a -> Frac b
frobenius_norm = |mx| sum(map_elementwise(mx, |x| x * x)) |> Num.to_frac |> Num.sqrt
