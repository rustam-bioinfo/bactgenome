def compute_stats(assembly):
    return {
        "total_contigs": len(assembly.contigs),
        "total_length": assembly.total_length(),
        "n50": assembly.n50(),
        "l50": assembly.l50(),
        "gc_content": round(assembly.gc_content() * 100, 2),
        "min_length": min((c.length() for c in assembly.contigs), default=0),
        "max_length": max((c.length() for c in assembly.contigs), default=0),
    }
