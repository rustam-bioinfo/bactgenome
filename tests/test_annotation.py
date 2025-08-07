from bactgenome.annotation import Annotation

def test_annotation_init():
    attributes = {"ID": "gene1", "Name": "geneA"}
    ann = Annotation(
        contig="scaffold1",
        source="Prodigal",
        feature_type="CDS",
        start=100,
        end=200,
        score=95.5,
        strand="+",
        phase=0,
        attributes=attributes,
    )
    assert ann.contig == "scaffold1"
    assert ann.source == "Prodigal"
    assert ann.feature_type == "CDS"
    assert ann.start == 100
    assert ann.end == 200
    assert ann.score == 95.5
    assert ann.strand == "+"
    assert ann.phase == 0
    assert ann.attributes == attributes

def test_annotation_length():
    ann = Annotation(
        contig="scaffold1",
        source="Prodigal",
        feature_type="CDS",
        start=100,
        end=200,
        score=None,
        strand="+",
        phase=None,
        attributes={},
    )
    assert ann.length == 101
