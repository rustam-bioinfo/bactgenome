from bactgenome.io import read_fasta, write_fasta
from bactgenome.qc import compute_stats

class Assembly:
    def __init__(self, contigs):
        self.contigs = contigs

    @classmethod
    def from_fasta(cls, file_path):
        return cls(read_fasta(file_path))

    def to_fasta(self, file_path):
        write_fasta(self.contigs, file_path)

    def total_length(self):
        return sum(c.length() for c in self.contigs)

    def n50(self):
        lengths = sorted((c.length() for c in self.contigs), reverse=True)
        total = sum(lengths)
        cum = 0
        for l in lengths:
            cum += l
            if cum >= total / 2:
                return l
        return 0

    def l50(self):
        lengths = sorted((c.length() for c in self.contigs), reverse=True)
        total = sum(lengths)
        cum, count = 0, 0
        for l in lengths:
            cum += l
            count += 1
            if cum >= total / 2:
                return count
        return 0

    def gc_content(self):
        total_len = self.total_length()
        total_gc = sum(c.gc_content() * c.length() for c in self.contigs)
        return total_gc / total_len if total_len else 0

    def contig_lengths(self):
        return [c.length() for c in self.contigs]

    def filter_by_length(self, min_len):
        self.contigs = [c for c in self.contigs if c.length() >= min_len]
        return self

    def sort_by_length(self, reverse=True):
        self.contigs.sort(key=lambda c: c.length(), reverse=reverse)
        return self

    def stats(self):
        return compute_stats(self)
