"""
Crosslinker definition functions.

Consider making these regexes in the future, or otherwise enhancing this module.
"""


def bs3(peptide1, peptide2):
    """
    BS3 can link any lysines, and also N-termini (not directly predicted here)
    """
    if "K" in peptide1 and "K" in peptide2:
        return True
    return False

