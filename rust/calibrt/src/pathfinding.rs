// --------------------------------------------------------------------------------
// Module 2: Optimal Ascending Path Identification
// --------------------------------------------------------------------------------

/// Minimum distance threshold for edge weight calculation.
/// Distances smaller than this are treated as degenerate to avoid numerical instability.
const DISTANCE_THRESHOLD: f64 = 1e-6;

/// The scratch buffers `find_optimal_path` reuses across calls, grouped so the
/// function itself doesn't take one parameter per buffer. Each one is cleared
/// and refilled, so repeated calls at a stable problem size allocate nothing.
#[derive(Default)]
pub(crate) struct PathfindingScratch {
    pub max_weights: Vec<f64>,
    pub prev_node_indices: Vec<Option<usize>>,
    /// Filled only when `ObserveOpts::dp_nodes` is set.
    pub considered: Vec<(usize, f64)>,
    /// The DP chain and the two greedy tails Pass 2 attaches, drained into the
    /// caller's `out_path` at the end of the call.
    pub path: Vec<crate::Point>,
    pub prefix: Vec<crate::Point>,
    pub suffix: Vec<crate::Point>,
}

/// Finds the highest-weight path through the nodes that satisfies the monotonic constraint.
///
/// This implements a dynamic programming solution on a directed acyclic graph (DAG) where:
/// - Nodes are sorted by (x, y) to ensure topological order
/// - Edges exist only between nodes where both x and y increase (monotonic constraint)
/// - Edge weights favor high-confidence nodes that are geometrically close
///
/// Writes the assembled path (DP chain plus any greedily-attached prefix/suffix,
/// see Pass 2 below) into `out_path`, clearing it first. Returns the
/// index range within that path the DP itself chose (`path[..range.start]`
/// and `path[range.end..]` are the greedy tails), and the DP recurrence's
/// objective value at the chosen end node (covers only `path[range]`).
pub(crate) fn find_optimal_path<O: crate::FitObserver>(
    nodes: &mut [crate::grid::Node],
    lookback: usize,
    scratch: &mut PathfindingScratch,
    out_path: &mut Vec<crate::Point>,
    obs: &mut O,
    opts: crate::ObserveOpts,
) -> (std::ops::Range<usize>, f64) {
    let PathfindingScratch {
        max_weights,
        prev_node_indices,
        considered,
        path,
        prefix,
        suffix,
    } = scratch;
    out_path.clear();
    path.clear();
    prefix.clear();
    suffix.clear();
    if nodes.is_empty() {
        return (0..0, 0.0);
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
        if opts.dp_nodes {
            considered.clear();
        }

        let start = i.saturating_sub(lookback);
        for j in start..i {
            // Only create edges where both dimensions increase (monotonic constraint)
            if nodes[i].center.library > nodes[j].center.library
                && nodes[i].center.observed > nodes[j].center.observed
            {
                let dx = nodes[i].center.library - nodes[j].center.library;
                let dy = nodes[i].center.observed - nodes[j].center.observed;
                let dist = (dx * dx + dy * dy).sqrt();

                if dist > DISTANCE_THRESHOLD {
                    // Edge weight formula: sqrt(weight_i) * sqrt(weight_j) / distance
                    // - Geometric mean of weights: Prefers high-confidence nodes but doesn't
                    //   annihilate edges to sparse-but-real cells (sqrt compresses the scale)
                    // - Division by distance: Penalizes long jumps, encouraging smooth curves
                    let edge_weight =
                        (nodes[i].center.weight.sqrt() * nodes[j].center.weight.sqrt()) / dist;
                    if opts.dp_nodes {
                        considered.push((j, edge_weight));
                    }
                    let new_weight = max_weights[j] + edge_weight;

                    if new_weight > max_weights[i] {
                        max_weights[i] = new_weight;
                        prev_node_indices[i] = Some(j);
                    }
                }
            }
        }

        if opts.dp_nodes {
            obs.on_event(crate::FitEvent::DpNode {
                i,
                node: &nodes[i],
                chose: prev_node_indices[i],
                acc_weight: max_weights[i],
                considered,
            });
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
    let mut current_idx_opt = Some(end_of_path_idx);
    while let Some(current_idx) = current_idx_opt {
        path.push(nodes[current_idx].center);
        current_idx_opt = prev_node_indices[current_idx];
    }
    path.reverse();

    if path.is_empty() {
        return (0..0, max_path_weight);
    }
    let dp_len = path.len();

    // Pass 2: Greedily extend the path beyond the DP endpoints.
    // The DP optimizes total weight and may skip sparse-but-valid regions at
    // the edges. We extend by walking through remaining non-suppressed nodes
    // that satisfy monotonicity, picking the nearest one at each step.

    // Extend backward: find nodes before the path start that satisfy monotonicity.
    // Nodes are sorted by (library, observed), so candidates are before the path's
    // first node in the sorted order.
    let first = path[0];
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
    if let Some(end_idx) = last_sorted_idx {
        let mut cursor = last;
        for node in nodes.iter().skip(end_idx + 1) {
            let candidate = node.center;
            if candidate.weight > 0.0
                && candidate.library > cursor.library
                && candidate.observed > cursor.observed
            {
                suffix.push(candidate);
                cursor = candidate;
            }
        }
    }

    // The DP's own span sits after `prefix` and before `suffix` in the
    // assembled path below.
    let dp_range = prefix.len()..(prefix.len() + dp_len);

    // Assemble into the caller-owned buffer: prefix + DP path + suffix.
    out_path.append(prefix);
    out_path.append(path);
    out_path.append(suffix);
    (dp_range, max_path_weight)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Node;
    use crate::{
        ObserveOpts,
        Point,
    };

    /// A `Node` at diagonal position `i` (library = observed = i + 0.5) with
    /// the given weight. `find_optimal_path` never reads `suppressed` or the
    /// `sum_*` accumulators, so this bypasses `Grid`/`suppress_nonmax`
    /// entirely — deliberately, since `suppress_nonmax` has its own floor
    /// (nodes below weight 1.0 with no competing row/column entry are
    /// dropped before ever reaching the DP) that would confound a fixture
    /// aimed at the DP's own accept/decline behavior.
    fn node_at(library: f64, observed: f64, weight: f64) -> Node {
        Node {
            center: Point {
                library,
                observed,
                weight,
            },
            suppressed: false,
            sum_wx: 0.0,
            sum_wy: 0.0,
            sum_w: 0.0,
        }
    }

    fn diag_node(i: usize, weight: f64) -> Node {
        let v = i as f64 + 0.5;
        node_at(v, v, weight)
    }

    /// Weights `[0.5, 10, 12, 11, 0.4]` on a diagonal: node 0 is a valid
    /// monotonic predecessor of node 1 but too weak to be worth routing
    /// through (chaining through it beats neither node 1's own weight of 10
    /// nor, transitively, anything reachable from node 1 onward), so the DP
    /// restarts fresh at node 1. Node 0 is still weight > 0 and monotonic, so
    /// Pass 2's backward walk attaches it as a prefix. Nodes 1..4 chain
    /// together under the DP (chaining through node 4 despite its own low
    /// weight of 0.4 still beats stopping at node 3, since the edge weight
    /// added is positive), so the DP's span is `1..5` of the assembled path.
    #[test]
    fn pass_two_attaches_a_dp_declined_leading_node_as_prefix() {
        let mut nodes: Vec<Node> = [0.5, 10.0, 12.0, 11.0, 0.4]
            .into_iter()
            .enumerate()
            .map(|(i, w)| diag_node(i, w))
            .collect();
        let mut scratch = PathfindingScratch::default();
        let mut path = Vec::new();

        let (dp_range, dp_weight) = find_optimal_path(
            &mut nodes,
            5,
            &mut scratch,
            &mut path,
            &mut (),
            ObserveOpts::NONE,
        );

        assert_eq!(
            path.len(),
            5,
            "all five nodes are monotonic and weight > 0, so all five appear \
             in the assembled path"
        );
        assert_eq!(
            dp_range,
            1..5,
            "the DP declines node 0, choosing only nodes 1..4"
        );
        assert!(dp_weight > 0.0);
    }

    /// The suffix mirror of the test above, which `dp_range` needs because
    /// widening `end` to `out_path.len()` swallows the whole greedy suffix while
    /// leaving a prefix-only test green.
    ///
    /// A trailing node the DP *includes* can never lower `max_weights`, so
    /// declining one takes a node with no admissible in-window predecessor: with
    /// `lookback == 1` the dip at `(3.5, 0.1)` fails the monotonic edge back to
    /// the core chain, and the stray at `(4.5, 4.5)` can only look back one rank
    /// — at the dip — so it accumulates the dip's 0.4 rather than the chain's,
    /// and the DP's best path ends at the chain. Pass 2's forward walk re-checks
    /// monotonicity against the DP's chosen *endpoint* instead, skips the dip
    /// (its observed RT is below the chain) and grafts the stray on as a suffix.
    #[test]
    fn pass_two_attaches_a_dp_declined_trailing_node_as_suffix() {
        let mut nodes = vec![
            diag_node(0, 10.0),
            diag_node(1, 12.0),
            diag_node(2, 11.0),
            node_at(3.5, 0.1, 0.4), // the dip
            node_at(4.5, 4.5, 0.5), // the stray
        ];
        let mut scratch = PathfindingScratch::default();
        let mut path = Vec::new();

        let (dp_range, _) = find_optimal_path(
            &mut nodes,
            1,
            &mut scratch,
            &mut path,
            &mut (),
            ObserveOpts::NONE,
        );

        assert_eq!(
            path.len(),
            4,
            "the three chain nodes plus the stray; the dip is not monotonic \
             against the chain and appears nowhere: {path:?}"
        );
        assert_eq!(
            dp_range,
            0..3,
            "the DP chooses only the core chain, so the stray is outside its \
             span: {path:?}"
        );
        assert!(
            dp_range.end < path.len(),
            "a `dp_range` reaching the end of the path would report the greedy \
             suffix as DP-chosen"
        );
        assert_eq!(
            path[dp_range.end].library, 4.5,
            "the one node past the DP's span must be the greedily attached stray"
        );
    }
}
