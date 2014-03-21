"""
Lightweight, simplistic memoization functions.
"""


def memoize_single(f):
    """ Memoization, works only on functions with one argument. """
    class MemoDict(dict):
        __slots__ = ()

        def __missing__(self, key):
            self[key] = ret = f(key)
            return ret

    return MemoDict().__getitem__


def memoize_args(f):
    """ Memoization supporting *args """
    class MemoDict(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            ret = self[key] = f(*key)
            return ret

    return MemoDict().__getitem__
