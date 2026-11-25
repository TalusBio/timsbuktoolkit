from scripts.gen_psims_origin_type import descendants, parse_terms


def test_decoy_spectrum_includes_transitive_subtypes():
    ontology = """[Term]
id: MS:1003192
name: decoy spectrum

[Term]
id: MS:9990002
name: test nested decoy subtype
is_a: MS:9990001 ! test decoy subtype

[Term]
id: MS:9990001
name: test decoy subtype
is_a: MS:1003192 ! decoy spectrum

[Term]
id: MS:1003073
name: observed spectrum
"""
    _, parents = parse_terms(ontology)

    assert descendants([1003192], parents) == {1003192, 9990001, 9990002}
