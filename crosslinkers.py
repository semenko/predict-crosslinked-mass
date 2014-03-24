"""
Crosslinker definition functions.

Each function should follow the format of:

Input:
    *peptides: 1 to N input peptides

Output:
    Dictionary of:
        can_link = bool True/False
        mass_shift = float (only if can_link is True)
"""


def bs3(*peptides):
    """
    BS3 can link any lysines, and also N-termini (not directly predicted here)
    """
    if len(peptides) == 1:
        # Mono-link
        if "K" in peptides[0]:
            return {'can_link': True, 'mass_shift': 156.07864}
    elif len(peptides) == 2:
        if "K" in peptides[0] and "K" in peptides[1]:
            return {'can_link': True, 'mass_shift': 138.06808}
    return {'can_link': False}
