# lau-topological-data-analysis

Topological data analysis (TDA) — extracting shape from data via persistent homology, simplicial complexes, and the Mapper algorithm.

## Features

- **Simplicial Complexes**: Vietoris-Rips, Čech, Alpha complexes, and Delaunay triangulation basics
- **Persistent Homology**: Filtration computation, barcode extraction, persistence diagrams
- **Betti Numbers**: Computed across filtrations at any scale
- **Distances**: Bottleneck and Wasserstein distances between persistence diagrams
- **Persistence Landscapes**: Construction, integration, and L^p distances
- **Mapper Algorithm**: Cluster-based simplification of high-dimensional data
- **Nerve Theorem**: Cover-based nerve construction and verification
- **Statistical TDA**: Bootstrap for persistence diagrams, confidence sets
- **Agent Behavior Analysis**: Topological features of agent state trajectories

## Usage

```rust
use lau_topological_data_analysis::*;
use nalgebra::DVector;

// Compute persistent homology
let points = vec![
    DVector::from_vec(vec![0.0, 0.0]),
    DVector::from_vec(vec![1.0, 0.0]),
    DVector::from_vec(vec![0.0, 1.0]),
];
let diagram = compute_persistent_homology(&points, 1);

// Get Betti numbers
let betti = betti_numbers(&diagram, 2.0);

// Compute bottleneck distance
let diagram2 = compute_persistent_homology(&points, 1);
let dist = bottleneck_distance(&diagram, &diagram2);

// Persistence landscape
let landscape = PersistenceLandscape::from_diagram(&diagram);
```

## License

MIT
