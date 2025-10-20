// --------------------------------------------------------------------------------
// Module 2: Optimal Ascending Path Identification
// --------------------------------------------------------------------------------

/// Finds the highest-weight path through the nodes that satisfies the monotonic constraint.
pub(crate) fn find_optimal_path(nodes: &mut [crate::grid::Node]) -> Vec<crate::Point> {
    if nodes.is_empty() {
        return Vec::new();
    }

    // Sort nodes primarily by x, then by y to process them in order for DAG pathfinding.
    nodes.sort_by(|a, b| {
        a.center
            .x
            .partial_cmp(&b.center.x)
            .unwrap()
            .then_with(|| a.center.y.partial_cmp(&b.center.y).unwrap())
    });

    let n = nodes.len();
    let mut max_weights = vec![0.0; n];
    let mut prev_node_indices = vec![None; n];

    for i in 0..n {
        max_weights[i] = nodes[i].center.weight; // Path can start at any node

        for j in 0..i {
            // 2.1 & 2.2: Check for monotonic edge and calculate weight
            if nodes[i].center.x > nodes[j].center.x && nodes[i].center.y > nodes[j].center.y {
                let dx = nodes[i].center.x - nodes[j].center.x;
                let dy = nodes[i].center.y - nodes[j].center.y;
                let dist = (dx * dx + dy * dy).sqrt();

                if dist > 1e-6 {
                    // Avoid division by zero
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

    for (i, max_w) in max_weights.into_iter().enumerate() {
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
