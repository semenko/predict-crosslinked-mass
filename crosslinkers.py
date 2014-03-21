"""
Crosslinker definition functions.

Each function should follow the format of:

Input:
    *peptides: 1-N input peptides

Output:
    Tuple of (bool, float), where:
        bool = can crosslink?
        float = mass shift (if linked, otherwise 0)
"""


def bs3(*peptides):
    """
    BS3 can link any lysines, and also N-termini (not directly predicted here)
    """
    if len(peptides) == 1:
        # Mono-link
        if "K" in peptides[0]:
            return (True, 156.07864)
    elif len(peptides) == 2:
        if "K" in peptides[0] and "K" in peptides[1]:
            return (True, 138.06808)
    return (False, 0)
