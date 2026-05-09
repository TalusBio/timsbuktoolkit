from bench._shuffle import generate_shuffled_entrapment, shuffle_keeping_c_terminal


def test_shuffle_keeps_c_terminal():
    """C-terminal residue stays put; interior is rearranged."""
    rng_seed = 42
    p = "PEPTIDEK"
    out = shuffle_keeping_c_terminal(p, seed=rng_seed)
    assert out[-1] == "K"
    # Same residues, possibly reordered
    assert sorted(out) == sorted(p)


def test_shuffle_short_peptide_returns_input():
    """A 2-AA peptide has no interior to shuffle; returns as-is."""
    out = shuffle_keeping_c_terminal("MK", seed=1)
    assert out == "MK"


def test_generate_shuffled_basic():
    targets = {"PEPTIDEK", "ANOTHERR"}
    pairs = generate_shuffled_entrapment(targets, r=1, seed=7)
    # One shuffle per target → 2 pairs
    assert len(pairs) == 2
    by_target = {t: s for (t, s) in pairs}
    assert "PEPTIDEK" in by_target
    assert "ANOTHERR" in by_target
    # Each shuffle preserves C-term + AA composition
    for t, s in pairs:
        assert s[-1] == t[-1]
        assert sorted(s) == sorted(t)
        # Shuffle is unique vs target peptides AND its own previously-generated shuffles
        assert s != t  # for non-trivially-permutable peptides
        assert s not in targets


def test_generate_shuffled_drops_target_when_no_unique_shuffle():
    """Target with fixed-point interior is dropped when r shuffles can't be found."""
    # 'AAAAK' interior 'AAAA' only permutes to itself → no unique shuffle exists.
    # 'AAK' interior 'A' has only one permutation.
    targets = {"AAAAK", "PEPTIDEK"}
    pairs = generate_shuffled_entrapment(targets, r=1, seed=42)
    matched_targets = {t for (t, _) in pairs}
    assert "AAAAK" not in matched_targets
    assert "PEPTIDEK" in matched_targets


def test_generate_shuffled_r_greater_than_one():
    targets = {"PEPTIDEK"}
    pairs = generate_shuffled_entrapment(targets, r=3, seed=11)
    # 3 shuffles for the one target
    assert len(pairs) == 3
    shuffles = [s for (_, s) in pairs]
    assert len(set(shuffles)) == 3  # all distinct
    for s in shuffles:
        assert s[-1] == "K"
        assert sorted(s) == sorted("PEPTIDEK")
        assert s != "PEPTIDEK"


def test_generate_shuffled_deterministic_with_seed():
    targets = {"PEPTIDEK", "ANOTHERR", "GREATSEQK"}
    a = generate_shuffled_entrapment(targets, r=2, seed=99)
    b = generate_shuffled_entrapment(targets, r=2, seed=99)
    assert a == b
