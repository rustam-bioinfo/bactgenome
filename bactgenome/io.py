from typing import List, Dict, Any
from .contig import Contig
from .annotation import Annotation

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

def _parse_attributes(attributes_str: str) -> Dict[str, Any]:
    """Parses the GFF attributes string (9th column)."""
    attributes = {}
    for part in attributes_str.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            attributes[key] = value
    return attributes

def read_gff(file_path: str) -> List[Annotation]:
    """
    Reads a GFF3 file and returns a list of Annotation objects.

    Args:
        file_path (str): The path to the GFF3 file.

    Returns:
        list[Annotation]: A list of Annotation objects.
    """
    annotations: List[Annotation] = []
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("##") or line.startswith("#!"):
                continue
            if line.strip() == "" or line.startswith("#"):
                continue

            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            contig = parts[0]
            source = parts[1]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])

            score_str = parts[5]
            score = float(score_str) if score_str != '.' else None

            strand = parts[6]

            phase_str = parts[7]
            phase = int(phase_str) if phase_str != '.' else None

            attributes_str = parts[8]
            attributes = _parse_attributes(attributes_str)

            annotations.append(
                Annotation(
                    contig=contig,
                    source=source,
                    feature_type=feature_type,
                    start=start,
                    end=end,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes=attributes,
                )
            )
    return annotations
