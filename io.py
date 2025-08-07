from typing import List
from bactgenome.contig import Contig

def read_fasta(file_path: str) -> List[Contig]:
    """
    Reads a FASTA file and returns a list of Contig objects.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        list[Contig]: A list of Contig objects.
    """
    contigs: List[Contig] = []
    with open(file_path, "r") as f:
        name, desc, seq_lines = None, "", []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    sequence = "".join(seq_lines)
                    contigs.append(Contig(name, sequence, desc))
                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                seq_lines = []
            else:
                seq_lines.append(line)
        if name:
            sequence = "".join(seq_lines)
            contigs.append(Contig(name, sequence, desc))
    return contigs

def write_fasta(contigs: List[Contig], file_path: str) -> None:
    """
    Writes a list of Contig objects to a FASTA file.

    Args:
        contigs (list[Contig]): A list of Contig objects to write.
        file_path (str): The path to the output FASTA file.
    """
    with open(file_path, "w") as f:
        for contig in contigs:
            header = f">{contig.name} {contig.description}\n".strip()
            f.write(header + "\n")
            for i in range(0, len(contig.sequence), 80):
                f.write(contig.sequence[i:i+80] + "\n")
