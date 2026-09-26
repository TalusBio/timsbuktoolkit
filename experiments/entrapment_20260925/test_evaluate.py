from experiments.entrapment_20260925.evaluate import evaluate


def test_selects_most_targets_under_empirical_limit_and_reports_fixed_q():
    rows = [
        {
            "q_value": "0.005",
            "paired_fdp": "0.008",
            "n_targets": "100",
            "n_entrapments": "1",
        },
        {
            "q_value": "0.009",
            "paired_fdp": "0.012",
            "n_targets": "120",
            "n_entrapments": "2",
        },
        {
            "q_value": "0.011",
            "paired_fdp": "0.009",
            "n_targets": "130",
            "n_entrapments": "3",
        },
    ]

    result = evaluate(rows)

    assert result["best_at_empirical_fdp"]["n_targets"] == 130
    assert result["at_reported_q"]["n_targets"] == 120
    assert result["at_reported_q"]["paired_fdp"] == 0.012
