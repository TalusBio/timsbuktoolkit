// --------------------------------------------------------------------------------
// Module 2: Optimal Ascending Path Identification
// --------------------------------------------------------------------------------

/// Minimum distance threshold for edge weight calculation.
/// Distances smaller than this are treated as degenerate to avoid numerical instability.
const DISTANCE_THRESHOLD: f64 = 1e-6;

/// Finds the highest-weight path through the nodes that satisfies the monotonic constraint.
///
/// This implements a dynamic programming solution on a directed acyclic graph (DAG) where:
/// - Nodes are sorted by (x, y) to ensure topological order
/// - Edges exist only between nodes where both x and y increase (monotonic constraint)
/// - Edge weights favor high-confidence nodes that are geometrically close
pub(crate) fn find_optimal_path(
    nodes: &mut [crate::grid::Node],
    lookback: usize,
    max_weights: &mut Vec<f64>,
    prev_node_indices: &mut Vec<Option<usize>>,
) -> Vec<crate::Point> {
    if nodes.is_empty() {
        return Vec::new();
    }

    // Sort nodes primarily by x, then by y to process them in order for DAG pathfinding.
    // This ensures we can use dynamic programming: when processing node i, all potential
    // predecessors (with smaller x,y) have already been processed.
    nodes.sort_by(|a, b| {
        a.center
            .x
            .partial_cmp(&b.center.x)
            .unwrap()
            .then_with(|| a.center.y.partial_cmp(&b.center.y).unwrap())
    });

    let n = nodes.len();
    max_weights.clear();
    max_weights.resize(n, 0.0);
    prev_node_indices.clear();
    prev_node_indices.resize(n, None);

    for i in 0..n {
        max_weights[i] = nodes[i].center.weight; // Path can start at any node

        let start = if i > lookback { i - lookback } else { 0 };
        for j in start..i {
            // Only create edges where both dimensions increase (monotonic constraint)
            if nodes[i].center.x > nodes[j].center.x && nodes[i].center.y > nodes[j].center.y {
                let dx = nodes[i].center.x - nodes[j].center.x;
                let dy = nodes[i].center.y - nodes[j].center.y;
                let dist = (dx * dx + dy * dy).sqrt();

                if dist > DISTANCE_THRESHOLD {
                    // Edge weight formula: (weight_i * weight_j) / distance
                    // - Product of weights: Prioritizes paths through high-confidence nodes
                    // - Division by distance: Penalizes long jumps, encouraging smooth curves
                    // This balances data fidelity (high weights) with geometric smoothness (short edges)
                    let edge_weight = (nodes[i].center.weight * nodes[j].center.weight) / dist;
                    let new_weight = max_weights[j] + edge_weight;

                    if new_weight > max_weights[i] {
                        max_weights[i] = new_weight;
                        prev_node_indices[i] = Some(j);
                    }
                }
            }
        }
    }

    // 2.3 Path Finding: Find the path with the maximum weight sum
    let mut max_path_weight = 0.0;
    let mut end_of_path_idx = 0;

    for (i, &max_w) in max_weights.iter().enumerate() {
        if max_w > max_path_weight {
            max_path_weight = max_w;
            end_of_path_idx = i;
        }
    }

    // Reconstruct the path
    let mut path = Vec::new();
    let mut current_idx_opt = Some(end_of_path_idx);
    while let Some(current_idx) = current_idx_opt {
        path.push(nodes[current_idx].center);
        current_idx_opt = prev_node_indices[current_idx];
    }
    path.reverse();

    path
}
