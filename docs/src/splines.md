# Splines

Isogeometric analysis uses B-splines and NURBS as the discrete approximation basis. A one-dimensional spline is built from a knot vector and a polynomial degree. Surfaces and solids are products of these one-dimensional functions. Control points set the shape. NURBS add a positive weight at each control point, which is what allows exact circles and other conic sections.

The patch type is a `NURBSMesh`. After Bézier extraction, see [Bezier extraction](@ref), assembly follows the [Ferrite documentation](https://ferrite-fem.github.io/Ferrite.jl/stable/).

## B-splines

A **knot vector** of degree $p$ with $n$ basis functions is a sequence of length $n+p+1$,

```math
\Xi = [\xi_1, \ldots, \xi_{n+p+1}].
```

Each interval $\xi_i < \xi_{i+1}$ is one element. Repeating a knot value raises its multiplicity $m$. At an interior knot of multiplicity $m$, with $1 \leq m \leq p+1$, a spline of degree $p$ is $C^{p-m}$. Derivatives through order $p-m$ match from the left and right. Higher derivatives may jump. A simple knot ($m=1$) gives $C^{p-1}$ continuity across that element boundary.

The B-spline basis is defined by the Cox–de Boor recursion. For $p=0$,

```math
\hat{N}_{A,0}(\xi) =
\begin{cases}
1 & \text{if } \xi_A \leq \xi < \xi_{A+1}, \\
0 & \text{otherwise},
\end{cases}
```

The functions sum to one, $\sum_A \hat{N}_{A,p}(\xi) = 1$, with the last knot interval closed on the right. For $p \geq 1$,

```math
\hat{N}_{A,p}(\xi) = \frac{\xi - \xi_A}{\xi_{A+p} - \xi_A} \, \hat{N}_{A,p-1}(\xi) + \frac{\xi_{A+p+1} - \xi}{\xi_{A+p+1} - \xi_{A+1}} \, \hat{N}_{A+1,p-1}(\xi).
```


A B-spline surface is a product of one-dimensional bases. With control points $\boldsymbol{X}_A$ and global index $A$ for the pair $(i,j)$,

```math
\boldsymbol{S}(\xi,\eta) = \sum_{A=1}^{N} \boldsymbol{X}_A \, N_A(\xi, \eta),
\qquad
N_A(\xi, \eta) = \hat{N}_{i}^{(\xi)}(\xi) \, \hat{N}_{j}^{(\eta)}(\eta).
```

The factors $\hat{N}_{i}^{(\xi)}$ and $\hat{N}_{j}^{(\eta)}$ use knot vectors $\Xi^{(\xi)}$ and $\Xi^{(\eta)}$ and degrees $p_\xi$, $p_\eta$.

### Open knot vectors

A knot vector is open for degree $p$ when the first and last values each have multiplicity $p+1$,

```math
\xi_1 = \cdots = \xi_{p+1},
\qquad
\xi_{n+1} = \cdots = \xi_{n+p+1}.
```

This is the standard choice in CAD and IGA. The first and last basis functions equal 1 at the ends, so the curve or surface passes through the first and last rows of control points. The active parameter range is $\xi \in [\xi_{p+1}, \xi_{n+1}]$, between the first and last distinct knots.

## NURBS

A **non-uniform rational B-spline** (NURBS) uses the same knot vector and control points, together with weights $w_A > 0$. The weight function is

```math
W(\xi) = \sum_{A=1}^{n} w_A \, \hat{N}_{A,p}(\xi),
```

and the rational basis is

```math
R_{A,p}(\xi) = \frac{w_A \, \hat{N}_{A,p}(\xi)}{W(\xi)}.
```

These functions sum to one, $\sum_A R_{A,p}(\xi) = 1$. If all weights are equal, $R_{A,p} = \hat{N}_{A,p}$.

A NURBS surface uses the same product of one-dimensional bases, with $R_A$ in place of $N_A$. With weights $w_{ij}$ on the control points,

```math
W(\xi,\eta) = \sum_{i=1}^{n_\xi} \sum_{j=1}^{n_\eta} w_{ij} \, \hat{N}_{i}^{(\xi)}(\xi) \, \hat{N}_{j}^{(\eta)}(\eta),
```

```math
R_A(\xi,\eta) = R_{ij}(\xi,\eta) = \frac{w_{ij} \, \hat{N}_{i}^{(\xi)}(\xi) \, \hat{N}_{j}^{(\eta)}(\eta)}{W(\xi,\eta)},
```

```math
\boldsymbol{S}(\xi,\eta) = \sum_{A=1}^{N} \boldsymbol{X}_A \, R_A(\xi, \eta).
```


## Refinement

The geometry is unchanged by the three operations below. They act on one parametric direction `dir` and rewrite copies of the knot vectors, control points, and weights. A new `NURBSMesh` is then built from those arrays.

**h-refinement** inserts a knot (Boehm). `knotinsertion!(kv, orders, cp, w, ξ; dir)` inserts the value `ξ` once in direction `dir`. The number of elements in that direction increases by one if `ξ` was not already a knot. Continuity at a newly inserted simple knot is $C^{p-1}$. Inserting a value that is already present raises the multiplicity and lowers the continuity to $C^{p-m}$.

**p-refinement** raises the polynomial degree by one (Piegl and Tiller, *The NURBS Book*, algorithm A5.9). `orders = orderelevation!(kv, orders, cp, w; dir)` returns the updated order tuple. End knots remain open (multiplicity $p+1$ at the new degree). Interior multiplicities increase by one as well, so $C^{p-m}$ is preserved.

**k-refinement** raises the degree first and then inserts new knots (Hughes, Cottrell, and Bazilevs, 2005). `orders = smoothnesselevation!(kv, orders, cp, w, new_knots; dir)` calls `orderelevation!` once and then `knotinsertion!` for each entry of `new_knots`. Those new knots have multiplicity one at the elevated degree, so continuity there is $C^{p}$. 

```julia
mesh = generate_nurbs_patch(:hypercube, (1, 1), (2, 2); cornerpos=(0.0, 0.0), size=(2.0, 3.0))

kv = (copy(mesh.knot_vectors[1]), copy(mesh.knot_vectors[2]))
cp = copy(mesh.control_points)
w  = copy(mesh.weights)
orders = mesh.orders

# h-refinement, insert ξ = 0 in the first parametric direction
knotinsertion!(kv, orders, cp, w, 0.0; dir=1)

# p-refinement, raise degree by one in that direction
orders = orderelevation!(kv, orders, cp, w; dir=1)

mesh = NURBSMesh(kv, orders, cp, w)
grid = BezierGrid(mesh)
```

Copy the arrays before refining. The knot vectors and control points are rewritten in place. Build a new `NURBSMesh` afterwards so the element connectivity is up to date.
