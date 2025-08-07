from bactgenome.io import read_gff

def test_read_gff():
    gff_file = "tests/data/prokka_example.gff"
    annotations = read_gff(gff_file)

    # There are 11 feature lines in the example file
    assert len(annotations) == 11

    # Test the first annotation in detail
    first_ann = annotations[0]
    assert first_ann.contig == "Scaffold_1_cov_18.34"
    assert first_ann.source == "Prodigal:002006"
    assert first_ann.feature_type == "CDS"
    assert first_ann.start == 10
    assert first_ann.end == 156
    assert first_ann.score is None
    assert first_ann.strand == "+"
    assert first_ann.phase == 0
    assert first_ann.attributes["ID"] == "NKBHOGAP_00001"
    assert first_ann.attributes["product"] == "Oxaloacetate decarboxylase beta chain"

def test_read_gff_empty_file(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "empty.gff"
    p.write_text("")
    annotations = read_gff(p)
    assert len(annotations) == 0

def test_read_gff_header_only(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "header.gff"
    p.write_text("##gff-version 3\n##sequence-region Scaffold_1_cov_18.34 1 476160")
    annotations = read_gff(p)
    assert len(annotations) == 0
