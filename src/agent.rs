//! Agent behavior analysis: topological features of agent state trajectories.

use crate::persistence::{compute_persistent_homology, PersistenceDiagram, PersistencePair, betti_numbers};
use crate::landscape::PersistenceLandscape;
use nalgebra::DVector;
use serde::{Serialize, Deserialize};

/// A topological feature extracted from agent behavior.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopologicalFeature {
    pub dim: usize,
    pub birth: f64,
    pub death: f64,
    pub persistence: f64,
    pub interpretation: String,
}

impl TopologicalFeature {
    pub fn from_pair(pair: &PersistencePair) -> Self {
        let interpretation = interpret_feature(pair);
        TopologicalFeature {
            dim: pair.dim,
            birth: pair.birth,
            death: pair.death,
            persistence: pair.persistence(),
            interpretation,
        }
    }
}

fn interpret_feature(pair: &PersistencePair) -> String {
    match pair.dim {
        0 => {
            if pair.is_essential() {
                "Connected component (essential)".to_string()
            } else {
                format!("Merged component at distance {:.3}, lived {:.3}", pair.death, pair.persistence())
            }
        }
        1 => format!("Loop/hole detected, radius range [{:.3}, {:.3}]", pair.birth, pair.death),
        2 => format!("Void/cavity detected, radius range [{:.3}, {:.3}]", pair.birth, pair.death),
        _ => format!("H{} feature, radius range [{:.3}, {:.3}]", pair.dim, pair.birth, pair.death),
    }
}

/// Analysis result for agent trajectories.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryAnalysis {
    /// Number of trajectory points.
    pub n_points: usize,
    /// Persistence diagram.
    pub diagram: PersistenceDiagram,
    /// Betti numbers at various scales.
    pub betti_at_scales: Vec<(f64, Vec<usize>)>,
    /// Extracted topological features.
    pub features: Vec<TopologicalFeature>,
    /// Persistence landscape.
    pub landscape: PersistenceLandscape,
    /// Summary metrics.
    pub summary: TrajectorySummary,
}

/// Summary metrics for trajectory analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectorySummary {
    /// Number of persistent H0 features (connected components that live long).
    pub persistent_components: usize,
    /// Number of persistent H1 features (loops).
    pub persistent_loops: usize,
    /// Maximum persistence.
    pub max_persistence: f64,
    /// Total persistence.
    pub total_persistence: f64,
    /// Estimated intrinsic dimensionality.
    pub estimated_dimension: usize,
}

/// Analyzer for agent state trajectories.
pub struct AgentTrajectoryAnalyzer {
    max_dim: usize,
}

impl AgentTrajectoryAnalyzer {
    pub fn new(max_dim: usize) -> Self {
        AgentTrajectoryAnalyzer { max_dim }
    }

    /// Analyze a trajectory of agent states.
    pub fn analyze(&self, states: &[Vec<f64>]) -> TrajectoryAnalysis {
        let points: Vec<DVector<f64>> = states.iter()
            .map(|s| DVector::from_vec(s.clone()))
            .collect();
        let n = points.len();

        let diagram = compute_persistent_homology(&points, self.max_dim);

        // Compute Betti numbers at multiple scales
        let dists = crate::complex::pairwise_distances(&points);
        let max_dist = dists.iter().fold(0.0f64, |a, &b| a.max(b));
        let scales: Vec<f64> = (1..=5).map(|i| max_dist * i as f64 / 5.0).collect();
        let betti_at_scales: Vec<(f64, Vec<usize>)> = scales.iter()
            .map(|&eps| (eps, betti_numbers(&diagram, eps)))
            .collect();

        let features: Vec<TopologicalFeature> = diagram.pairs.iter()
            .map(|p| TopologicalFeature::from_pair(p))
            .collect();

        let landscape = PersistenceLandscape::from_diagram(&diagram);

        // Compute summary
        let persistent_threshold = max_dist * 0.1;
        let persistent_components = diagram.pairs_of_dim(0).iter()
            .filter(|p| p.persistence() > persistent_threshold || p.is_essential())
            .count();
        let persistent_loops = diagram.pairs_of_dim(1).iter()
            .filter(|p| p.persistence() > persistent_threshold)
            .count();
        let max_pers = diagram.max_persistence();
        let total_pers = diagram.total_persistence();

        // Estimate dimension from Betti numbers
        let max_betti_dim = betti_at_scales.iter()
            .map(|(_, bn)| bn.len())
            .max()
            .unwrap_or(1)
            .saturating_sub(1);

        TrajectoryAnalysis {
            n_points: n,
            diagram,
            betti_at_scales,
            features,
            landscape,
            summary: TrajectorySummary {
                persistent_components,
                persistent_loops,
                max_persistence: max_pers,
                total_persistence: total_pers,
                estimated_dimension: max_betti_dim.max(1),
            },
        }
    }

    /// Compare two trajectories by their topological signatures.
    pub fn compare(&self, states1: &[Vec<f64>], states2: &[Vec<f64>]) -> TrajectoryComparison {
        let points1: Vec<DVector<f64>> = states1.iter()
            .map(|s| DVector::from_vec(s.clone()))
            .collect();
        let points2: Vec<DVector<f64>> = states2.iter()
            .map(|s| DVector::from_vec(s.clone()))
            .collect();

        let dg1 = compute_persistent_homology(&points1, self.max_dim);
        let dg2 = compute_persistent_homology(&points2, self.max_dim);

        let bottleneck = crate::distance::bottleneck_distance(&dg1, &dg2);
        let wasserstein = crate::distance::wasserstein_distance(&dg1, &dg2, 2.0);

        let l1 = PersistenceLandscape::from_diagram(&dg1);
        let l2 = PersistenceLandscape::from_diagram(&dg2);
        let landscape_dist = crate::landscape::landscape_distance(&l1, &l2, 2.0);

        TrajectoryComparison {
            bottleneck_distance: bottleneck,
            wasserstein_distance: wasserstein,
            landscape_distance: landscape_dist,
            similar: bottleneck < 0.5,
        }
    }
}

/// Result of comparing two trajectories.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryComparison {
    pub bottleneck_distance: f64,
    pub wasserstein_distance: f64,
    pub landscape_distance: f64,
    pub similar: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analyze_simple_trajectory() {
        let states: Vec<Vec<f64>> = (0..4).map(|i| {
            vec![i as f64, 0.0]
        }).collect();
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let analysis = analyzer.analyze(&states);
        assert_eq!(analysis.n_points, 4);
        assert!(!analysis.features.is_empty());
        assert!(analysis.summary.persistent_components >= 1);
    }

    #[test]
    fn test_analyze_circular_trajectory() {
        // Points on a circle (small)
        let states: Vec<Vec<f64>> = (0..4).map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / 4.0;
            vec![theta.cos(), theta.sin()]
        }).collect();
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let analysis = analyzer.analyze(&states);
        assert_eq!(analysis.n_points, 4);
        // Should detect a loop (H1 feature)
        assert!(analysis.features.iter().any(|f| f.dim == 1));
    }

    #[test]
    fn test_compare_same_trajectory() {
        let states: Vec<Vec<f64>> = (0..4).map(|i| vec![i as f64, 0.0]).collect();
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let comp = analyzer.compare(&states, &states);
        assert!(comp.bottleneck_distance < 1e-10);
        assert!(comp.similar);
    }

    #[test]
    fn test_compare_different_trajectories() {
        let line: Vec<Vec<f64>> = (0..4).map(|i| vec![i as f64, 0.0]).collect();
        let circle: Vec<Vec<f64>> = (0..4).map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / 4.0;
            vec![theta.cos(), theta.sin()]
        }).collect();
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let comp = analyzer.compare(&line, &circle);
        assert!(comp.bottleneck_distance > 0.0);
    }

    #[test]
    fn test_topological_feature_interpretation() {
        let pair = PersistencePair::new(0, 0.0, 1.0);
        let feature = TopologicalFeature::from_pair(&pair);
        assert_eq!(feature.dim, 0);
        assert!((feature.persistence - 1.0).abs() < 1e-10);
        assert!(!feature.interpretation.is_empty());
    }

    #[test]
    fn test_trajectory_analysis_serialization() {
        let states: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64, 0.0]).collect();
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let analysis = analyzer.analyze(&states);
        let json = serde_json::to_string(&analysis).unwrap();
        assert!(json.len() > 10);
    }

    #[test]
    fn test_analyze_single_state() {
        let states = vec![vec![1.0, 2.0]];
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let analysis = analyzer.analyze(&states);
        assert_eq!(analysis.n_points, 1);
        assert!(analysis.summary.persistent_components >= 1);
    }

    #[test]
    fn test_analyze_two_states() {
        let states = vec![vec![0.0, 0.0], vec![1.0, 0.0]];
        let analyzer = AgentTrajectoryAnalyzer::new(1);
        let analysis = analyzer.analyze(&states);
        assert_eq!(analysis.n_points, 2);
    }

    #[test]
    fn test_trajectory_comparison_serialization() {
        let comp = TrajectoryComparison {
            bottleneck_distance: 0.5,
            wasserstein_distance: 0.3,
            landscape_distance: 0.2,
            similar: true,
        };
        let json = serde_json::to_string(&comp).unwrap();
        let c2: TrajectoryComparison = serde_json::from_str(&json).unwrap();
        assert!(c2.similar);
    }
}
