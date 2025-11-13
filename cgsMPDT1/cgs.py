# Bismillahirrahmanirrahim

"""
About Module
------------
This module is a Python implementation of the CGS Theory with CGS spaces equipped with CGS-MPD-T1.
"""


# -------------------------------------------------------------------
# PYTHON MODULES IN USE
import numpy as np
from random import choices
from math import floor


# -------------------------------------------------------------------
# MAX-PLUS SCRIPT
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


# -------------------------------------------------------------------
# CGS-MPD-T1 SCRIPT

class CGSMPDT1:
    def __init__(self, lenPAC, lenSAC, n, d):
        self.lenPAC = lenPAC
        self.lenSAC = lenSAC
        self.PAC = [f'p{k}' for k in range(1, lenPAC + 1)]
        self.SAC = [f'p{k}' for k in range(1, lenSAC + 1)]
        self.n = n
        self.d = d
        d1, d2, d3, d4, d5 = d
        self.Delta = oprod([d1, d2, d3]) |ot| (d4|op|d5)
        self.__inv_arg = "ERROR: Invalid arguments!"

    def __isinPAC(self, m):
        try:
            if isinstance(m, str) \
                    and m[0] == 'p' \
                    and isinstance(int(m[1]), int):
                return True
            else:
                return False
        except:
            return False
        
    def __isinSAC(self, m):
        try:
            if isinstance(m, str) \
                    and m[0] == 's' \
                    and isinstance(int(m[1]), int):
                return True
            else:
                return False
        except:
            return False

    # Ordering of PAC:
    def omega_p(self, m):
        if self.__isinPAC(m):
            return int(m[1])
        else:
            print(self.__inv_arg)
            raise ValueError
        
    # Inversed of Omega_p:
    def omega_pInv(self, i):
        try:
            if isinstance(i, int) and 0 < i <= self.lenPAC:
                return f'p{i}'
            else:
                print(self.__inv_arg)
                raise ValueError
        except:
            print(self.__inv_arg)
            raise ValueError
        
    # Ordering of SAC:
    def omega_s(self, m):
        if self.__isinSAC(m):
            return int(m[1])
        else:
            print(self.__inv_arg)
            raise ValueError
        
    # Inversed of Omega_s:
    def omega_sInv(self, i):
        try:
            if isinstance(i, int) and 0 < i <= self.lenSAC:
                return f's{i}'
            else:
                print(self.__inv_arg)
                raise ValueError
        except:
            print(self.__inv_arg)
            raise ValueError
        
    # Utilization function:
    def rho(self, m):
        n1, n2 = self.n
        lenPAC, lenSAC = self.lenPAC, self.lenSAC
        if self.__isinPAC(m):
            r = (n1 *n2) % lenPAC
            q = (n1 *n2) // lenPAC
            i = self.omega_p(m)
        
        elif self.__isinSAC(m):
            r = n2 % lenSAC
            q = n2 // lenSAC
            i = self.omega_s(m)

        else:
            print(self.__inv_arg)
            raise ValueError
        
        return q + 1 if 0 < i <= r else q
    
    # Primary Grid:
    def primGRID(self):
        return [
            (m, k) for m in self.PAC
            for k in range(1, self.rho(m) + 1)
            ]
    
    # Far Right Stack:
    def secGRID(self):
        return [
            (m, k) for m in self.SAC
            for k in range(1, self.rho(m) + 1)
            ]
    
    # Entire Rectangle:
    def entGRID(self):
        calEp = self.primGRID()
        calEs = self.secGRID()
        return calEp + calEs
    
    # Lexicographic ordering:
    def lexi(self, b):
        calEp = self.primGRID()
        calEs = self.secGRID()
        lenPAC, lenSAC = self.lenPAC, self.lenSAC
        n1, n2 = self.n
        
        if b in calEp:
            m, k = b
            t = self.omega_p(m) + (k - 1) *lenPAC
            return t + floor((t - 1) / n1)
        elif b in calEs:
            m, k = b
            return (n1 + 1) *(self.omega_s(m) + (k - 1) *lenSAC)
        else:
            print(self.__inv_arg)
            raise ValueError
        
    # Inversed Lexicographic ordering:
    def antiLexi(self, a):
        n1, n2 = self.n
        lenPAC, lenSAC = self.lenPAC, self.lenSAC
        try:
            if isinstance(a, int) and a <= (n1 + 1) *n2:
                pass
            else:
                print(self.__inv_arg)
                raise ValueError
        except:
            print(self.__inv_arg)
            raise ValueError
        
        if a % (n1 + 1) > 0:
            k = floor((a - floor(a / (n1 + 1)) -1) / lenPAC) + 1
            m = self.omega_pInv(
                a - floor(a / (n1 + 1)) - (k - 1) *lenPAC
            )
        elif a % (n1 + 1) == 0:
            k = floor(( a/(n1 + 1) -1 ) / lenSAC) + 1
            m = self.omega_sInv(
                int(a / (n1 + 1)) - (k - 1) *lenSAC
            )
        else:
            print(self.__inv_arg)
            raise ValueError
        
        return (m, k)
    
    # Auxilliary function for generating CGS-DM-T1
    def __maxplus_aux(self, c_num):
        """
        Note
        ----
        The scripts in this function follows from Defintion 17 of the original manuscript.
        """
        n1, n2 = self.n
        d1, d2, d3, d4, _ = self.d
        dim = (n1 + 1) *n2
        R = [[E for k in range(dim)] for k in range(dim)]
        
        # Axiom C1
        if c_num == 1:
            for p in range(dim):
                for q in range(dim):
                    i, j = p + 1, q + 1
                    if i == j + 1 and i % (n1 + 1) != 1:
                        R[p][q] = d1 |ot| 1
        
        # Axiom C2
        elif c_num == 2:
            for p in range(dim):
                for q in range(dim):
                    i, j = p + 1, q + 1
                    if i == j + (n1 + 1):
                        R[p][q] = oprod([d1, d2, d3]) |ot| 1

        # Axiom C3
        elif c_num == 3:
            for p in range(dim):
                for q in range(dim):
                    i, j = p + 1, q + 1
                    if i == j + 2 *(n1 + 1):
                        R[p][q] = oprod([d1, d2, d3, d4])

        # Axiom C4
        elif c_num == 4:
            for p in range(dim):
                for q in range(dim):
                    i, j = p + 1, q + 1
                    mi, ki = self.antiLexi(i)
                    mj, kj = self.antiLexi(j)

                    if mi == mj and ki == kj + 1:
                        R[p][q] = oprod([d1, d2, d3, d4])
        return R
    
    # CGS-DM-T1
    def cgsDMT1(self, triangularity= False):
        R_list = [self.__maxplus_aux(k) for k in range(1, 5)]
        # Axiom C5:
        R = osum(R_list)

        # Lower Matrix-Triangularity test for CGS-DM-T1:
        if triangularity == False:
            return R
        elif triangularity == True:
            truth = list()
            for i in range(len(R)):
                for j in range(len(R[i])):
                    if j > i and R[i][j] == E:
                        truth.append(True)
            return R, all(truth)
        else:
            print(self.__inv_arg)
            raise ValueError
    
    # CGS State Closure Operator:
    def cgsSCO(self, u= None):
        n1, n2 = self.n
        dim = (n1 + 1) *n2
        if u is None:
            u = [E] *dim
            u[0] = 0
        R = self.cgsDMT1()
        # Follows from Definition 11 in the paper:
        return (R|ast|dim) |ot| u
    
    # CGS Functional:
    def cgsF(self, u= None):
        # Follows from Definition 18:
        Delta = self.Delta
        # Follows from Theorem 9:
        return Delta |ot| osum(self.cgsSCO(u))

    # CGS Primary Functional
    def cgsPF(self, u= None):
        Delta = self.Delta
        x = self.cgsSCO(u)[-2]
        return Delta |ot| x
    
    # CGS Secondary Functional:
    def cgsSF(self, u= None):
        n1, _ = self.n
        Delta = self.Delta
        x, y = self.cgsSCO()[n1], osum(self.cgsSCO(u))
        return (Delta |ot| y) |ot_inv| x
    
    # CGS Primary-Secondary Functional
    def cgs_allF(self, u= None):
        n1, _ = self.n
        Delta = self.Delta
        sco = self.cgsSCO(u)
        v = osum(sco)
        w = sco[-2]
        x = sco[n1]
        return Delta |ot| v, Delta |ot| w, (Delta |ot| v) |ot_inv| x