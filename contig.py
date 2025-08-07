class Contig:
    """Represents a single contig in a genome assembly."""

    def __init__(self, name: str, sequence: str, description: str = "") -> None:
        """
        Initializes a Contig object.

        Args:
            name (str): The name of the contig.
            sequence (str): The DNA sequence of the contig.
            description (str, optional): A description of the contig. Defaults to "".
        """
        self.name: str = name
        self.sequence: str = sequence.upper()
        self.description: str = description

    def length(self) -> int:
        """
        Calculates the length of the contig.

        Returns:
            int: The length of the contig's sequence.
        """
        return len(self.sequence)

    def gc_content(self) -> float:
        """
        Calculates the GC content of the contig.

        Returns:
            float: The GC content as a fraction.
        """
        g = self.sequence.count("G")
        c = self.sequence.count("C")
        return (g + c) / len(self.sequence) if self.sequence else 0

    def reverse_complement(self) -> str:
        """
        Computes the reverse complement of the contig's sequence.

        Returns:
            str: The reverse complemented sequence.
        """
        complement = str.maketrans("ATCG", "TAGC")
        return self.sequence.translate(complement)[::-1]
