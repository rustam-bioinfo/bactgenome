from typing import Dict, Any
from bactgenome.assembly import Assembly

def compute_stats(assembly: Assembly) -> Dict[str, Any]:
    """
    Computes a dictionary of statistics for a given assembly.

    Args:
        assembly (Assembly): The assembly to compute statistics for.

    Returns:
        dict: A dictionary of assembly statistics.
    """
    return {
        "total_contigs": len(assembly.contigs),
        "total_length": assembly.total_length(),
        "n50": assembly.n50(),
        "l50": assembly.l50(),
        "gc_content": round(assembly.gc_content() * 100, 2),
        "min_length": min((c.length() for c in assembly.contigs), default=0),
        "max_length": max((c.length() for c in assembly.contigs), default=0),
    }
