from typing import Dict, Any, Optional

class Annotation:
    """Represents a single feature from a GFF file."""

    def __init__(
        self,
        contig: str,
        source: str,
        feature_type: str,
        start: int,
        end: int,
        score: Optional[float],
        strand: str,
        phase: Optional[int],
        attributes: Dict[str, Any],
    ) -> None:
        """
        Initializes an Annotation object.

        Args:
            contig (str): The name of the sequence where the feature is located.
            source (str): The source of the feature (e.g., a program like Prodigal).
            feature_type (str): The type of the feature (e.g., 'CDS', 'gene').
            start (int): The starting position of the feature (1-based).
            end (int): The ending position of the feature (1-based).
            score (Optional[float]): The score of the feature.
            strand (str): The strand of the feature ('+' or '-').
            phase (Optional[int]): The phase for CDS features (0, 1, or 2).
            attributes (Dict[str, Any]): A dictionary of attributes from the 9th column of the GFF file.
        """
        self.contig = contig
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __repr__(self) -> str:
        return (
            f"Annotation(contig='{self.contig}', type='{self.feature_type}', "
            f"start={self.start}, end={self.end}, strand='{self.strand}')"
        )

    @property
    def length(self) -> int:
        """Calculates the length of the feature."""
        return self.end - self.start + 1
