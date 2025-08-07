import pytest
from bactgenome.assembly import Assembly
from bactgenome.contig import Contig
from bactgenome.io import read_gff
from bactgenome.annotation import Annotation

@pytest.fixture
def assembly_with_annotations():
    # A contig long enough to contain all features from the test GFF
    contigs = [Contig("Scaffold_1_cov_18.34", "A" * 9000)]
    gff_file = "tests/data/prokka_example.gff"
    annotations = read_gff(gff_file)
    return Assembly(contigs, annotations=annotations)

def test_add_annotations_from_gff():
    # Create a dummy assembly
    contigs = [Contig("Scaffold_1_cov_18.34", "ATGC")]
    assembly = Assembly(contigs)

    # Check that it has no annotations initially
    assert len(assembly.annotations) == 0

    # Add annotations from the example GFF file
    gff_file = "tests/data/prokka_example.gff"
    assembly.add_annotations_from_gff(gff_file)

    # Check that the annotations have been added
    assert len(assembly.annotations) == 11

    # Check a detail from the first annotation to be sure
    first_ann = assembly.annotations[0]
    assert first_ann.attributes["ID"] == "NKBHOGAP_00001"

def test_init_with_annotations():
    # Create a dummy assembly with annotations in the constructor
    contigs = [Contig("Scaffold_1_cov_18.34", "ATGC")]
    gff_file = "tests/data/prokka_example.gff"
    from bactgenome.io import read_gff
    annotations = read_gff(gff_file)

    assembly = Assembly(contigs, annotations=annotations)

    # Check that the annotations have been added
    assert len(assembly.annotations) == 11
    first_ann = assembly.annotations[0]
    assert first_ann.attributes["ID"] == "NKBHOGAP_00001"

def test_get_annotation_by_id(assembly_with_annotations):
    # Test getting an existing annotation
    ann = assembly_with_annotations.get_annotation_by_id("NKBHOGAP_00003")
    assert ann is not None
    assert ann.attributes["product"] == "Serine endoprotease DegS"

    # Test getting a non-existent annotation
    ann = assembly_with_annotations.get_annotation_by_id("non_existent_id")
    assert ann is None

def test_get_sequence_for_annotation(assembly_with_annotations):
    # Test a feature on the '+' strand
    ann_pos = assembly_with_annotations.get_annotation_by_id("NKBHOGAP_00001") # pos 10-156 on + strand
    seq = assembly_with_annotations.get_sequence_for_annotation(ann_pos)
    assert len(seq) == 147
    assert seq == "A" * 147

    # Test a feature on the '-' strand
    ann_neg = assembly_with_annotations.get_annotation_by_id("NKBHOGAP_00003") # pos 952-2010 on - strand
    seq = assembly_with_annotations.get_sequence_for_annotation(ann_neg)
    assert len(seq) == 1059
    assert seq == "T" * 1059 # Reverse complement of 'A's is 'T's

def test_get_sequence_for_annotation_not_found(assembly_with_annotations):
    # Create an annotation on a contig that doesn't exist
    ann = Annotation("wrong_scaffold", "test", "CDS", 1, 10, None, "+", None, {})
    with pytest.raises(ValueError):
        assembly_with_annotations.get_sequence_for_annotation(ann)
