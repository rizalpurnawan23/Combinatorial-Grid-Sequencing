# Bismillahirrahmanirrahim
# 2026-05-14

# --------------------------------------------------------------
# PYTHON MODULES IN USE
import numpy as np
from random import choices
from math import floor


# --------------------------------------------------------------
# PART 1: MAX-PLUS SCRIPT
class Infix:
    def __init__(self, func):
        self.func = func
    def __ror__(self, other):
        return Infix(
            lambda x,
            self= self,
            other= other: self.func(other, x)
            )
    def __or__(self, other):
        return self.func(other)
    
# E = 'E'     # Oplus identity
E = float("-inf")

# Identity matrix generator:
def id(n:int):
    """
    Parameters
    ----------
    n   : an integer representing the dimension of the identity
          max-plus matrix
    """
    idm = [[E for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                idm[i][j] = 0
    return idm

# Generate a matrix randomly:
def maxplus_random(m:int, n= None):
    """
    Parameters
    ----------
    m   : an integer representing the dimension of the columns
    n   : an integer representing the dimension of the rows,
          or a None

    Trivia: the max-plus semifield entries in this matrix are
    only integers in between -100 and 100 inclusive or the
    max-plus zero element (E).
    """
    if n is None:
        n = m
    candidates = list(range(-100, 101)) + [E] *23
    M = list()
    for i in range(m):
        M.append(choices(candidates, k= n))
    return M

# Auxilliary functions:
def __is_MPmatrix(a):
    try:
        return isinstance(a, list) \
            and all(isinstance(row, list) for row in a) \
                and all(
                    isinstance(x, int)
                    or isinstance(x, float)
                    or x == E
                    for row in a for x in row
                    )
    except:
        return False
    
def __is_MPsemimodule(a):
    try:
        return isinstance(a, list) \
                and all(
                    isinstance(x, int)
                    or isinstance(x, float)
                    or x == E
                    for x in a)
    except:
        return False
    
def __is_maxplus(a):
    if isinstance(a, int) or isinstance(a, float) or a == E:
        return True
    else:
        return False

# Matrix Transpose Function:
def tp(A):
    """
    Parameters
    ----------
    A   : a max-plus matrix
    """
    return [[x[k] for x in A] for k in range(len(A[0]))]

# Oplus Function
@Infix
def op(a, b):
    """
    Parameters
    ----------
    a   : a max-plus semifield element, max-plus semimodule
          element, or a max-plus matrix
    b   : a max-plus semifield element, max-plus semimodule
          element, or a max-plus matrix

    Trivia:
    > If both a and b are max-plus semimodule elements, they
      must be of the same dimension.
    > If both a and b are max-plus matrices, then the row
      dimension of a must equal to the column dimension of b.
    """
    # For max-plus matrices:
    if __is_MPmatrix(a) and __is_MPmatrix(b):
        x, y = np.array(a), np.array(b)
        c = list()
        for x, y in zip(a, b):
            d = list()
            for p, q in zip(x, y):
                d.append(max(p, q))
            c.append(d)
        return c
    
    # For max-plus semimodule elements:
    elif __is_MPsemimodule(a) and __is_MPsemimodule(b):
        x, y = np.array(a), np.array(b)
        c = list()
        for x, y in zip(a, b):
            c.append(max(x, y))
        return c
    
    # For semifield elements:
    elif __is_maxplus(a) and __is_maxplus(b):
        return max(a, b)
    
    else:
        print("ERROR: Invalid arguments!")
        raise ValueError
    
# Oplus Summation
def osum(li):
    """
    Parameters
    ----------
    li  : a list of either max-plus semifield elements, or
          max-plus semimodule elements of the same dimension,
          or max-plus matrices of the same dimension
    """
    # For max-plus matrices:
    if __is_MPmatrix(li[0]):
        mat = list()
        m, n = len(li[0]), len(li[0][0])
        for i in range(m):
            row = list()
            for j in range(n):
                row.append(max([x[i][j] for x in li]))
            mat.append(row)
        return mat
    
    # For max-plus semimodule elements:
    elif __is_MPsemimodule(li[0]):
        u = list()
        dim = len(li[0])
        for k in range(dim):
            u.append(max([x[k] for x in li]))
        return u
    
    # For max-plus elements:
    elif __is_maxplus(li[0]):
        return max(li)
    
    else:
        print("ERROR: Invalid arguments!")
        raise ValueError


# Otimes Function        
@Infix
def ot(a, b):
    """
    Parameters
    ----------
    a   : a max-plus semifield element, max-plus semimodule
          element, or a max-plus matrix
    b   : a max-plus semifield element, max-plus semimodule
          element, or a max-plus matrix

    Trivia:
    > If both a and b are max-plus semimodule elements, they
      must be of the same dimension.
    > If a is a max-plus matrix and b is a max-plus semimodule
      element, then the row dimension of a must equal to the
      dimension of b.
    > If both a and b are max-plus matrices, then the row
      dimension of a must equal to the column dimension of b.
    """
    # For matrix-matrix multiplication:
    if __is_MPmatrix(a) and __is_MPmatrix(b):
        c = list()
        for h in tp(b):
            d = list()
            for g in a:
                d.append(max([x + y for x, y in zip(g, h)]))
            c.append(d)
        return tp(c)
    
    # For matrix-semimodule multiplication:
    elif __is_MPmatrix(a) and __is_MPsemimodule(b):        
        c = list()
        for g in a:
            c.append(max([x + y for x, y in zip(g, b)]))
        return c

    # For max-plus semifield:
    elif __is_maxplus(a) and __is_maxplus(b):
        return a + b
    
    else:
        print("ERROR: Invalid arguments!")
        raise ValueError
    
@Infix
def ot_inv(a, b):
    if __is_maxplus(a) and __is_maxplus(b):
        pass
    else:
        print("ERROR: Invalid arguments!")
        raise ValueError
    return a - b
    
# Otimes Product
def oprod(li):
    """
    Parameters
    ----------
    li  : a list of either max-plus semifield elements, or
          max-plus semimodule elements of the same dimension,
          or max-plus matrices of the same dimension
    """
    # For max-plus matrices:
    if __is_MPmatrix(li[0]):
        x = li[0]
        for y in li[1:]:
            x = x |ot| y
        return x
    
    # For max-plus semifield elements:
    elif __is_maxplus(li[0]):
        return sum(li)
    
# Exponentiation
@Infix
def exp(a, q):
    """
    Parameters
    ----------
    a   : A max-plus matrix, or a max-plus semimodule element, or
          a max-plus semifield element.
    q   : a nonnegative integer.
    """
    # For square matrices:
    if __is_MPmatrix(a):
        i = 0
        x = id(len(a))
        while i < q:
            x = x |ot| a
            i += 1
        return x
  
    elif __is_maxplus(a):
        return q *a

    else:
        print("ERROR: Invalid arguments!")
        raise ValueError
        
# Truncated Kleene-star:
@Infix
def ast(A, n):
    """
    Parameters
    ----------
    A   : a square max-plus matrix or a max-plus semifield elemnt
    n   : a nonnegative integer
    """
    bundle = [A|exp|k for k in range(n)]
    return osum(bundle)


# --------------------------------------------------------------
# PART 2: CGS-MPD IMPLEMENTATION

class CGSSimpleStructure:
    def __init__(self, lenPAC, lenSAC, n, d):
        pass