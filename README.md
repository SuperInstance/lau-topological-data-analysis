# lau-topological-data-analysis

> Topological Data Analysis (TDA) in Rust — simplicial complexes, persistent homology, persistence landscapes, Mapper algorithm, diagram distances, bootstrap confidence, and agent trajectory analysis.

## What This Does

`lau-topological-data-analysis` extracts the **shape** of data using tools from algebraic topology. It builds simplicial complexes from point clouds, computes persistent homology to identify topological features (connected components, loops, voids), and provides statistical methods like bootstrap confidence sets for assessing feature significance. A dedicated agent module applies these tools to analyze behavioral trajectories.

The library covers the full TDA pipeline:

1. **Complexes** — Vietoris-Rips, Čech, and Alpha complexes with Delaunay triangulation.
2. **Persistence** — Filtration, barcode computation via matrix reduction, Betti numbers.
3. **Distances** — Bottleneck and Wasserstein distances between persistence diagrams.
4. **Landscapes** — Persistence landscapes with integration and Lᵖ norms.
5. **Mapper** — Cluster-based topological simplification of high-dimensional data.
6. **Nerve** — Nerve theorem verification for covers.
7. **Statistics** — Bootstrap resampling and confidence sets for persistence diagrams.
8. **Agent** — Topological analysis of agent state trajectories.

---

## Key Idea

Data has **shape**. Persistent homology captures that shape at every scale simultaneously:

- At small scales, everything is separate (many components).
- At intermediate scales, clusters merge and loops appear.
- At large scales, everything connects into one blob.

The features that persist across a wide range of scales are the **signal**; the short-lived ones are **noise**. This library computes those persistence pairs, visualizes them as diagrams or landscapes, and provides rigorous statistical tools for comparing them.

---

## Install

Add to your `Cargo.toml`:

```toml
[dependencies]
lau-topological-data-analysis = "0.1"
```

Requires **Rust 2021 edition** or later. Dependencies: `nalgebra`, `serde`, `serde_json`, `itertools`, `rand`.

---

## Quick Start

```rust
use lau_topological_data_analysis::*;
use nalgebra::DVector;

// 1. Create a point cloud
let points: Vec<DVector<f64>> = vec![
    DVector::from_vec(vec![0.0, 0.0]),
    DVector::from_vec(vec![1.0, 0.0]),
    DVector::from_vec(vec![0.0, 1.0]),
    DVector::from_vec(vec![1.0, 1.0]),
];

// 2. Compute persistent homology
let diagram = persistence::compute_persistent_homology(&points, 2);

// 3. Extract features
for pair in &diagram.pairs {
    println!("H{}: born at {:.3}, died at {:.3} (persistence: {:.3})",
        pair.dim, pair.birth, pair.death, pair.persistence());
}

// 4. Betti numbers at a scale
let betti = persistence::betti_numbers(&diagram, 1.0);
println!("Betti numbers at ε=1.0: {:?}", betti); // [β₀, β₁]

// 5. Persistence landscape
let landscape = landscape::PersistenceLandscape::from_diagram(&diagram);
println!("Landscape layers: {}", landscape.num_layers());
let value = landscape.evaluate(0, 0.5);
let integral = landscape.integrate(0);
```

---

## API Reference

### `complex` — Simplicial Complexes

#### Simplex

```rust
let s = Simplex::new(vec![0, 1, 2]);  // 2-simplex (triangle)
s.dim();           // 2
s.vertices();      // [0, 1, 2]
s.faces();         // [Simplex({0,1}), Simplex({0,2}), Simplex({1,2})]
s.contains_face(&Simplex::edge(0, 1));  // true
```

| Constructor | Description |
|-------------|-------------|
| `Simplex::new(vertices)` | From an iterator of vertex indices. |
| `Simplex::vertex(i)` | 0-simplex (point). |
| `Simplex::edge(i, j)` | 1-simplex (line segment). |

| Method | Description |
|--------|-------------|
| `dim() → usize` | Dimension (vertices - 1). |
| `vertices() → Vec<usize>` | Sorted vertex list. |
| `faces() → Vec<Simplex>` | Codimension-1 boundary simplices. |
| `contains_face(&Simplex) → bool` | Subface check. |

#### SimplicialComplex

A set of simplices closed under taking faces.

```rust
let sc = SimplicialComplex::from_simplices(vec![
    Simplex::new(vec![0, 1, 2])
]);
sc.n_simplices();  // 7 (1 triangle + 3 edges + 3 vertices)
sc.euler_characteristic();  // 1 (3 - 3 + 1)
```

| Method | Description |
|--------|-------------|
| `new()` | Empty complex. |
| `from_simplices(simplices)` | Build with automatic face closure. |
| `add_simplex(&simplex)` | Add a simplex and all its faces. |
| `simplices()` | Iterator over all simplices. |
| `simplices_of_dim(d)` | Filter by dimension. |
| `n_vertices()`, `n_simplices()` | Counts. |
| `contains(&simplex)` | Membership test. |
| `euler_characteristic() → i64` | χ = Σ (-1)^d × (number of d-simplices). |
| `boundary_matrix(dim) → DMatrix<i32>` | Z₂ boundary matrix for dimension dim. |

#### VietorisRips

```rust
let vr = VietorisRips::new(points, epsilon, max_dim);
let complex = vr.build();
vr.distance_matrix();  // Pairwise Euclidean distances
```

A simplex is included iff all pairwise distances between vertices ≤ epsilon.

#### CechComplex

```rust
let cech = CechComplex::new(points, epsilon, max_dim);
let complex = cech.build();
```

A simplex is included iff the intersection of balls of radius epsilon around all vertices is non-empty (diameter/2 ≤ epsilon condition).

#### AlphaComplex

```rust
let alpha = AlphaComplex::new(points, alpha_value);
let complex = alpha.build();
```

Delaunay-based complex filtered by circumradius. Full 2D Delaunay triangulation with empty circumcircle property; higher dimensions use distance-based approximation.

#### Utilities

```rust
let dists = pairwise_distances(&points);  // DMatrix<f64>
```

### `persistence` — Persistent Homology

#### PersistencePair

```rust
let p = PersistencePair::new(1, 0.5, 2.0);  // H₁ feature
p.persistence();  // 1.5
p.is_essential(); // false (finite death)
p.midpoint();     // 1.25
```

#### PersistenceDiagram

```rust
let dg = PersistenceDiagram::new(pairs);
dg.pairs_of_dim(0);       // H₀ features only
dg.len();                 // Total pairs
dg.max_persistence();     // Longest-lived feature
dg.total_persistence();   // Sum of all persistences
```

#### Filtration

```rust
let mut filt = Filtration::new();
filt.add(0.0, Simplex::vertex(0));
filt.add(1.5, Simplex::edge(0, 1));
filt.sort();
let diagram = filt.compute_persistence();
```

Or build a Vietoris-Rips filtration directly:

```rust
let filt = Filtration::vietoris_rips_filtration(&points, max_dim);
let diagram = filt.compute_persistence();
```

**Algorithm:** Standard matrix reduction over Z₂ (persistent homology via column operations).

#### Top-level Functions

| Function | Description |
|----------|-------------|
| `compute_persistent_homology(points, max_dim)` | One-shot VR filtration + persistence. |
| `betti_numbers(diagram, epsilon)` | Betti numbers at filtration value epsilon. |

### `distance` — Diagram Distances

#### Bottleneck Distance

```rust
let dist = bottleneck_distance(&dg1, &dg2);
```

The minimum ε such that all points in one diagram can be matched to points (or diagonal) in the other within L^∞ distance ε. Computed via binary search + maximum bipartite matching.

#### Wasserstein Distance

```rust
let dist = wasserstein_distance(&dg1, &dg2, 2.0);
```

Lᵖ-optimal matching between diagrams (Hungarian algorithm for small instances, greedy for larger).

### `landscape` — Persistence Landscapes

```rust
let landscape = PersistenceLandscape::from_diagram(&diagram);
landscape.num_layers();           // Number of λ_k functions
landscape.evaluate(k, t);         // λ_k(t)
landscape.integrate(k);           // ∫ λ_k(t) dt
landscape.norm_p(2.0);            // L² norm across all layers
```

A persistence landscape converts a diagram into a sequence of piecewise-linear functions. Each persistence pair becomes a "tent function" and layers are built by taking the k-th largest tent at each point.

**Landscape Distance:**

```rust
let dist = landscape_distance(&l1, &l2, 2.0);
```

### `mapper` — Mapper Algorithm

```rust
let config = MapperConfig {
    n_intervals: 10,
    overlap: 0.1,
    n_clusters: 3,
};
let mapper = Mapper::new(config);
let graph = mapper.run(&points, |p| p[0]);  // Filter by first coordinate

graph.n_nodes();
graph.n_edges();
graph.is_connected();
graph.n_components();
graph.adjacency();  // HashMap<usize, Vec<usize>>
```

**MapperConfig:**
- `n_intervals`: Number of intervals in the filter function domain.
- `overlap`: Fraction of overlap between intervals (0.0–1.0).
- `n_clusters`: Clusters per interval (k-means).

The Mapper algorithm:
1. Apply a filter function `f: ℝᵈ → ℝ` to all points.
2. Divide the filter range into overlapping intervals.
3. Cluster points within each interval (k-means).
4. Connect clusters that share points → Mapper graph.

### `nerve` — Nerve Theorem

```rust
let cover = Cover::new(vec![vec![0, 1], vec![1, 2], vec![0, 2]], 3);
cover.is_valid();   // true (covers all 3 points)
let nerve = cover.nerve();  // SimplicialComplex

let verification = verify_nerve_theorem(&cover);
// verification.cover_is_valid
// verification.contractibility_assumed
// verification.nerve

// Convert Mapper graph to cover
let cover = mapper_graph_to_cover(&graph, n_points);
```

The nerve of a cover is a simplicial complex where a simplex is included iff the corresponding sets have non-empty intersection. By the Nerve Theorem, if all sets and intersections are contractible, the nerve is homotopy equivalent to the union.

### `statistics` — Bootstrap & Confidence

```rust
// Bootstrap resampling
let result = statistics::bootstrap_persistence(&points, 50, 2, 0.95);
// result.diagrams: Vec<PersistenceDiagram>
// result.mean_bottleneck: average bottleneck distance to original
// result.std_bottleneck: std deviation
// result.confidence_band: mean + z × std

// Confidence set
let cs = statistics::confidence_set(&points, 2, 0.95, 100);
cs.center;    // Original diagram
cs.radius;    // Confidence band radius
cs.level;     // 0.95
cs.contains(&other_diagram);  // bottleneck ≤ radius?
```

Uses normal-approximation z-scores for confidence levels (1.96 for 95%, 2.576 for 99%).

### `agent` — Trajectory Analysis

```rust
let analyzer = AgentTrajectoryAnalyzer::new(2);  // max_dim = 2

// Analyze a single trajectory
let states: Vec<Vec<f64>> = vec![
    vec![0.0, 0.0], vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0],
];
let analysis = analyzer.analyze(&states);
analysis.n_points;
analysis.diagram;         // PersistenceDiagram
analysis.betti_at_scales; // Vec<(f64, Vec<usize>)>
analysis.features;        // Vec<TopologicalFeature> (with human-readable interpretations)
analysis.landscape;       // PersistenceLandscape
analysis.summary.persistent_components;
analysis.summary.persistent_loops;
analysis.summary.max_persistence;
analysis.summary.total_persistence;
analysis.summary.estimated_dimension;

// Compare two trajectories
let comp = analyzer.compare(&states1, &states2);
comp.bottleneck_distance;
comp.wasserstein_distance;
comp.landscape_distance;
comp.similar;  // true if bottleneck < 0.5
```

**TopologicalFeature** provides human-readable interpretations:
- H₀: "Connected component (essential)" or "Merged component at distance X"
- H₁: "Loop/hole detected, radius range [X, Y]"
- H₂: "Void/cavity detected, radius range [X, Y]"

---

## How It Works

### Full TDA Pipeline

```
Point Cloud
    │
    ├── Build Simplicial Complex (VR / Čech / Alpha)
    │       │
    │       └── Filtration (increasing epsilon)
    │               │
    │               └── Matrix Reduction over Z₂
    │                       │
    │                       └── Persistence Diagram
    │                               │
    │                               ├── Betti Numbers
    │                               ├── Persistence Landscape
    │                               ├── Bottleneck / Wasserstein Distance
    │                               └── Bootstrap Confidence
    │
    └── Mapper Algorithm
            │
            └── Mapper Graph → Nerve
```

### Vietoris-Rips Complex

For n points with pairwise distances, include a k-simplex iff all (k+1 choose 2) pairwise distances ≤ ε. Built by clique enumeration: first add vertices, then edges within ε, then extend to higher simplices.

### Persistent Homology (Matrix Reduction)

Given a filtration (simplices sorted by birth time), the algorithm:

1. Build the boundary matrix over Z₂.
2. Reduce columns from left to right using column operations.
3. Each reduced column's lowest non-zero row gives a persistence pair (birth, death).
4. Unpaired columns are **essential** features (infinite persistence).

### Persistence Landscapes

Convert each persistence pair (b, d) into a tent function:
```
λ(t) = max(0, min(t - b, d - t))
```
Layer k takes the (k+1)-th largest tent value at each point. This gives a stable, functional representation suitable for statistical analysis.

### Mapper Algorithm

1. Choose a filter function f: ℝᵈ → ℝ (e.g., first coordinate, distance from centroid, PCA projection).
2. Divide the filter range into overlapping intervals.
3. Within each interval, cluster points using k-means.
4. Build a graph: nodes = clusters, edges = shared points between clusters in adjacent intervals.

### Bootstrap for Persistence

1. Compute original persistence diagram D₀.
2. Resample points with replacement → compute diagram Dᵢ.
3. Compute bottleneck distance d(D₀, Dᵢ).
4. Report mean ± z×std of distances as confidence band.

---

## The Math

**Euler Characteristic:**

$$\chi = \sum_{d=0}^{D} (-1)^d \cdot |\{d\text{-simplices}\}|$$

**Boundary Matrix (over Z₂):**

$$[\partial_d]_{ij} = \begin{cases} 1 & \text{if } \sigma_i^{(d-1)} \text{ is a face of } \sigma_j^{(d)} \\ 0 & \text{otherwise} \end{cases}$$

**Betti Numbers at Scale ε:**

$$\beta_d(\varepsilon) = |\{(b, d') \in \text{Dgm} : b \leq \varepsilon < d'\}| + |\{(b, \infty) : b \leq \varepsilon\}|$$

**Bottleneck Distance:**

$$d_B(D_1, D_2) = \inf_{\text{matchings } \gamma} \sup_{x \in D_1} \|x - \gamma(x)\|_\infty$$

Where unmatched points are projected to the diagonal.

**Wasserstein-p Distance:**

$$W_p(D_1, D_2) = \left(\inf_{\gamma} \sum_{x \in D_1} \|x - \gamma(x)\|_p^p\right)^{1/p}$$

**Persistence Landscape:**

$$\lambda_k(t) = \text{k-th largest } \Lambda_{(b,d)}(t) \text{ where } \Lambda_{(b,d)}(t) = \max(0, \min(t-b, d-t))$$

**Landscape Lᵖ Norm:**

$$\|\lambda\|_p = \left(\sum_k \int |\lambda_k(t)|^p dt\right)^{1/p}$$

**Vietoris-Rips Condition:**

A k-simplex {v₀, ..., vₖ} is included iff ∀i,j: d(vᵢ, vⱼ) ≤ ε.

**Čech Condition:**

A k-simplex is included iff the intersection of ε-balls around all vertices is non-empty (diameter/2 ≤ ε).

---

## Tests

77 unit tests across all modules:

- **complex (13 tests):** simplex construction, faces, complex from simplices, Euler characteristic, Vietoris-Rips at various epsilon, Čech complex, Alpha complex, boundary matrix, pairwise distances, face containment
- **persistence (12 tests):** persistence pair properties, essential features, diagram properties, filtration (1/2/3 points, collinear), Betti numbers, symmetric difference, serialization
- **distance (9 tests):** bottleneck (identical, empty, shifted, different size, symmetry), Wasserstein (identical, p=2), diagonal distance, permutation
- **landscape (10 tests):** construction, evaluate peak/outside, integrate, empty diagram, multiple features, landscape distance (identical), tent evaluation, piecewise integration, serialization
- **mapper (9 tests):** basic, empty, single point, connected graph, clusters, serialization, components, disconnected, cluster_points
- **nerve (9 tests):** valid/invalid cover, nerve (basic, disjoint, triple overlap), verify theorem, mapper-to-cover, empty/single set
- **statistics (6 tests):** bootstrap (basic, single point), confidence set, z-score, serialization
- **agent (9 tests):** analyze (simple, circular, single, two states), compare (same, different), feature interpretation, serialization

Run with:
```bash
cargo test
```

---

## License

MIT
