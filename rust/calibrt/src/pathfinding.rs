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
            .library
            .partial_cmp(&b.center.library)
            .unwrap()
            .then_with(|| a.center.observed.partial_cmp(&b.center.observed).unwrap())
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
            if nodes[i].center.library > nodes[j].center.library && nodes[i].center.observed > nodes[j].center.observed {
                let dx = nodes[i].center.library - nodes[j].center.library;
                let dy = nodes[i].center.observed - nodes[j].center.observed;
                let dist = (dx * dx + dy * dy).sqrt();

                if dist > DISTANCE_THRESHOLD {
                    // Edge weight formula: sqrt(weight_i) * sqrt(weight_j) / distance
                    // - Geometric mean of weights: Prefers high-confidence nodes but doesn't
                    //   annihilate edges to sparse-but-real cells (sqrt compresses the scale)
                    // - Division by distance: Penalizes long jumps, encouraging smooth curves
                    let edge_weight = (nodes[i].center.weight.sqrt() * nodes[j].center.weight.sqrt()) / dist;
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

    // Reconstruct the DP path
    let mut path = Vec::new();
    let mut current_idx_opt = Some(end_of_path_idx);
    while let Some(current_idx) = current_idx_opt {
        path.push(nodes[current_idx].center);
        current_idx_opt = prev_node_indices[current_idx];
    }
    path.reverse();

    if path.is_empty() {
        return path;
    }

    // Pass 2: Greedily extend the path beyond the DP endpoints.
    // The DP optimizes total weight and may skip sparse-but-valid regions at
    // the edges. We extend by walking through remaining non-suppressed nodes
    // that satisfy monotonicity, picking the nearest one at each step.

    // Extend backward: find nodes before the path start that satisfy monotonicity.
    // Nodes are sorted by (library, observed), so candidates are before the path's
    // first node in the sorted order.
    let first = path[0];
    let mut prefix = Vec::new();
    // Walk backward through sorted nodes, greedily picking the nearest monotonic predecessor
    let first_sorted_idx = nodes.iter().position(|n| {
        (n.center.library - first.library).abs() < 1e-9
            && (n.center.observed - first.observed).abs() < 1e-9
    });
    if let Some(start_idx) = first_sorted_idx {
        let mut cursor = first;
        for j in (0..start_idx).rev() {
            let candidate = nodes[j].center;
            if candidate.weight > 0.0
                && candidate.library < cursor.library
                && candidate.observed < cursor.observed
            {
                prefix.push(candidate);
                cursor = candidate;
            }
        }
        prefix.reverse();
    }

    // Extend forward: find nodes after the path end that satisfy monotonicity.
    let last = *path.last().unwrap();
    let last_sorted_idx = nodes.iter().rposition(|n| {
        (n.center.library - last.library).abs() < 1e-9
            && (n.center.observed - last.observed).abs() < 1e-9
    });
    let mut suffix = Vec::new();
    if let Some(end_idx) = last_sorted_idx {
        let mut cursor = last;
        for j in (end_idx + 1)..n {
            let candidate = nodes[j].center;
            if candidate.weight > 0.0
                && candidate.library > cursor.library
                && candidate.observed > cursor.observed
            {
                suffix.push(candidate);
                cursor = candidate;
            }
        }
    }

    // Assemble: prefix + DP path + suffix
    let mut full_path = prefix;
    full_path.append(&mut path);
    full_path.append(&mut suffix);
    full_path
}
