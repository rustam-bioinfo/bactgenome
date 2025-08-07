from __future__ import annotations
from typing import List, Dict, Any

from bactgenome.contig import Contig
from bactgenome.io import read_fasta, write_fasta
from bactgenome.qc import compute_stats

class Assembly:
    """Represents a collection of contigs that make up a genome assembly."""

    def __init__(self, contigs: List[Contig]) -> None:
        """
        Initializes an Assembly object.

        Args:
            contigs (list[Contig]): A list of Contig objects.
        """
        self.contigs: List[Contig] = contigs

    @classmethod
    def from_fasta(cls, file_path: str) -> Assembly:
        """
        Creates an Assembly object from a FASTA file.

        Args:
            file_path (str): The path to the FASTA file.

        Returns:
            Assembly: A new Assembly object.
        """
        return cls(read_fasta(file_path))

    def to_fasta(self, file_path: str) -> None:
        """
        Writes the assembly to a FASTA file.

        Args:
            file_path (str): The path to the output FASTA file.
        """
        write_fasta(self.contigs, file_path)

    def total_length(self) -> int:
        """
        Calculates the total length of all contigs in the assembly.

        Returns:
            int: The total length.
        """
        return sum(c.length() for c in self.contigs)

    def n50(self) -> int:
        """
        Calculates the N50 statistic for the assembly.

        Returns:
            int: The N50 value.
        """
        lengths = sorted((c.length() for c in self.contigs), reverse=True)
        total = sum(lengths)
        cum = 0
        for l in lengths:
            cum += l
            if cum >= total / 2:
                return l
        return 0

    def l50(self) -> int:
        """
        Calculates the L50 statistic for the assembly.

        Returns:
            int: The L50 value.
        """
        lengths = sorted((c.length() for c in self.contigs), reverse=True)
        total = sum(lengths)
        cum, count = 0, 0
        for l in lengths:
            cum += l
            count += 1
            if cum >= total / 2:
                return count
        return 0

    def gc_content(self) -> float:
        """
        Calculates the overall GC content of the assembly.

        Returns:
            float: The GC content as a fraction.
        """
        total_len = self.total_length()
        total_gc = sum(c.gc_content() * c.length() for c in self.contigs)
        return total_gc / total_len if total_len else 0

    def contig_lengths(self) -> List[int]:
        """
        Gets a list of all contig lengths.

        Returns:
            list[int]: A list of contig lengths.
        """
        return [c.length() for c in self.contigs]

    def filter_by_length(self, min_len: int) -> Assembly:
        """
        Filters contigs in the assembly by a minimum length.

        Args:
            min_len (int): The minimum length of contigs to keep.

        Returns:
            Assembly: The current Assembly object with filtered contigs.
        """
        self.contigs = [c for c in self.contigs if c.length() >= min_len]
        return self

    def sort_by_length(self, reverse: bool = True) -> Assembly:
        """
        Sorts the contigs in the assembly by length.

        Args:
            reverse (bool, optional): Whether to sort in descending order. Defaults to True.

        Returns:
            Assembly: The current Assembly object with sorted contigs.
        """
        self.contigs.sort(key=lambda c: c.length(), reverse=reverse)
        return self

    def stats(self) -> Dict[str, Any]:
        """
        Computes a dictionary of statistics for the assembly.

        Returns:
            dict: A dictionary of assembly statistics.
        """
        return compute_stats(self)
