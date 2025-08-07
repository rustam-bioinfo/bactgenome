from bactgenome.assembly import Assembly
from bactgenome.contig import Contig

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
