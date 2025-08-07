class Contig:
    def __init__(self, name, sequence, description=""):
        self.name = name
        self.sequence = sequence.upper()
        self.description = description

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        g = self.sequence.count("G")
        c = self.sequence.count("C")
        return (g + c) / len(self.sequence) if self.sequence else 0

    def reverse_complement(self):
        complement = str.maketrans("ATCG", "TAGC")
        return self.sequence.translate(complement)[::-1]
