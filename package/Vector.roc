## A module for mathematical vectors with various operations.
module [
    Vector,
    new,
    multiply,
    scale,
    add_scalar,
    subtract_scalar,
    count,
    element,
    to_list,
    add,
    subtract,
    concat,
    to_frac,
    approx_eq,
    sum,
    mean,
    map,
    map2,
    distance,
    divide_elementwise,
    multiply_elementwise,
    normalize,
    magnitude,
    l1_norm,
    l2_norm,
    linf_norm,
    zeros,
    ones,
]

Vector a := {
    length : U64,
    values : List (Num a),
}
    implements [Eq, Inspect]

## Create a new vector from a list of numbers.
new : List (Num a) -> Vector a
new = |l|
    @Vector(
        {
            length: List.len(l),
            values: l,
        },
    )

## Get the number of elements in the vector.
count : Vector a -> U64
count = |@Vector(v)| v.length

## Get the element at the specified index.
element : Vector a, U64 -> Result (Num a) [OutOfBounds]
element = |@Vector(v), index| List.get(v.values, index)

## Convert the vector to a list of its elements.
to_list : Vector a -> List (Num a)
to_list = |@Vector(v)| v.values

## Element-wise addition of two vectors.
add : Vector a, Vector a -> Result (Vector a) [VectorLengthMismatch]
add = |@Vector(vec1), @Vector(vec2)|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(vec1.values, vec2.values, |v1, v2| v1 + v2)
        |> new
        |> Ok

## Subtract two vectors.
subtract : Vector a, Vector a -> Result (Vector a) [VectorLengthMismatch]
subtract = |@Vector(vec1), @Vector(vec2)|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(vec1.values, vec2.values, |v1, v2| v1 - v2)
        |> new
        |> Ok

## Dot product of two vectors.
multiply : Vector a, Vector a -> Result (Num a) [VectorLengthMismatch]
multiply = |@Vector(vec1), @Vector(vec2)|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(vec1.values, vec2.values, |v1, v2| v1 * v2) |> List.sum |> Ok

## Scale the vector by a scalar.
scale : Vector a, Num a -> Vector a
scale = |@Vector(vec), n| List.map(vec.values, |v| Num.mul(v, n)) |> new

## Add a scalar to each element of the vector.
add_scalar : Vector a, Num a -> Vector a
add_scalar = |@Vector(vec), n| List.map(vec.values, |v| v + n) |> new

## Subtract a scalar from each element of the vector.
subtract_scalar : Vector a, Num a -> Vector a
subtract_scalar = |@Vector(vec), n| List.map(vec.values, |v| v - n) |> new

## Concatenate two vectors.
concat : Vector a, Vector a -> Vector a
concat = |@Vector(vec1), @Vector(vec2)| List.concat(vec1.values, vec2.values) |> new

## Convert a vector of numbers to a vector of fractions.
to_frac : Vector * -> Vector (FloatingPoint *)
to_frac = |@Vector(vec)| List.map(vec.values, Num.to_frac) |> new

## Approximate equality check between two vectors containing floating-point numbers.
approx_eq : Vector (FloatingPoint a), Vector (FloatingPoint a) -> Bool
approx_eq = |@Vector(vec1), @Vector(vec2)|
    List.map2(vec1.values, vec2.values, |v1, v2| Num.is_approx_eq(v1, v2, {})) |> List.all(|b| b)

## The sum of the vector elements.
sum : Vector a -> Num a
sum = |@Vector(vec)| List.sum(vec.values)

## The mean (average) of the vector elements.
mean : Vector * -> Frac *
mean = |@Vector(vec)| List.sum(vec.values) |> Num.to_frac |> Num.div(Num.to_frac(vec.length))

## Element-wise division of two vectors.
divide_elementwise : Vector (FloatingPoint a), Vector (FloatingPoint a) -> Result (Vector (FloatingPoint a)) [VectorLengthMismatch, DivByZero]
divide_elementwise = |@Vector(vec1), @Vector(vec2)|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(vec1.values, vec2.values, |v1, v2| Num.div_checked(v1, v2))
        |> List.map_try(|v| v)
        |> Result.map_ok(new)

## Element-wise multiplication of two vectors.
multiply_elementwise : Vector a, Vector a -> Result (Vector a) [VectorLengthMismatch]
multiply_elementwise = |@Vector(vec1), @Vector(vec2)|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(vec1.values, vec2.values, |v1, v2| v1 * v2) |> new |> Ok

## Element-wise map of a vector with a unary function.
map : Vector a, (Num a -> Num b) -> Vector b
map = |@Vector(vec), f| List.map(vec.values, f) |> new

## Element-wise map of two vectors with a binary function.
map2 : Vector a, Vector b, (Num a, Num b -> Num c) -> Result (Vector c) [VectorLengthMismatch]
map2 = |@Vector(vec1), @Vector(vec2), f|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(vec1.values, vec2.values, f) |> new |> Ok

## The Euclidean distance between two vectors.
distance : Vector a, Vector a -> Result (Frac b) [VectorLengthMismatch]
distance = |@Vector(vec1), @Vector(vec2)|
    if vec1.length != vec2.length then
        Err(VectorLengthMismatch)
    else
        List.map2(
            vec1.values,
            vec2.values,
            |v1, v2|
                diff = Num.sub(v1, v2)
                Num.mul(diff, diff),
        )
        |> List.sum
        |> Num.to_frac
        |> Num.sqrt
        |> Ok

## The L1 norm, or manhattan norm. Sum of the absolute values of the vector elements.
l1_norm : Vector a -> Num a
l1_norm = |@Vector(vec)| List.map(vec.values, |v| Num.abs(v)) |> List.sum

## The L2 norm, or Euclidean norm. Square root of the sum of the squares of the vector elements.
l2_norm : Vector a -> Frac b
l2_norm = |@Vector(vec)|
    List.map(vec.values, |v| Num.mul(v, v))
    |> List.sum
    |> Num.to_frac
    |> Num.sqrt

## The L-infinity norm, or maximum norm. The maximum absolute value among the vector elements.
linf_norm : Vector a -> Num a
linf_norm = |@Vector(vec)| List.map(vec.values, Num.abs) |> List.max |> Result.with_default(0)

## The magnitude of the vector, defined as its L2 norm. Also referred to as the Euclidean length.
magnitude : Vector a -> Frac b
magnitude = |vec| l2_norm(vec)

## Normalize the vector to have a magnitude of 1. Returns an error if the vector has zero magnitude.
normalize : Vector a -> Result (Vector (FloatingPoint b)) [DivByZero]
normalize = |@Vector(vec)|
    length = magnitude(@Vector(vec))
    if Num.is_approx_eq(length, Num.to_frac(0), {}) then
        Err(DivByZero)
    else
        List.map(vec.values, |v| Num.div(Num.to_frac(v), length))
        |> new
        |> Ok

## Create a vector of the given length filled with zeros.
zeros : U64 -> Vector a
zeros = |len| List.repeat(0, len) |> new

## Create a vector of the given length filled with ones.
ones : U64 -> Vector a
ones = |len| List.repeat(1, len) |> new
