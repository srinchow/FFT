"""
# @author Srinjoy Choudhury <srinchow@gmail.com>
# @date 27/4/2020

# The Fast Polynomial Multiplication Algorithm
# A implementation of FFT to multiply 2 polynomials 
# Input: A,B, two arrays of integers representing polynomials
#        their length is in O(n)
# Output: Coefficient representation of AB
# Complexity: O(n logn)


"""


from cmath import exp
from math import pi


class NthRootOfUnity:
    def __init__(self, n, k = 1):
        self.k = k
        self.n = n

    def __pow__(self, other):
        if type(other) is int:
            n = NthRootOfUnity(self.n, self.k * other)
            return n

    def __eq__(self, other):
        if other == 1:
            return abs(self.n) == abs(self.k)

    def __mul__(self, other):
        return exp(2*1j*pi*self.k/self.n)*other

    def __repr__(self):
        return str(self.n) + "-th root of unity to the " + str(self.k)

    @property
    def th(self):
        return abs(self.n // self.k)


class FPM:
    def __init__(self,A,B):
        self.A = A
        self.B = B
        self.size = 1<<(len(A)+len(B)-2).bit_length()
    
    def FFT(self,A,w):
        if w == 1:
            return [sum(A)]
        o2 = w**2
        Aeven = self.FFT(A[0::2], o2)
        Aodd = self.FFT(A[1::2], o2)
        C3 = [None]*w.th
        for i in range(w.th//2):
            C3[i] = Aeven[i] + w**i * Aodd[i]
            C3[i+w.th//2] = Aeven[i] - w**i * Aodd[i]
        return C3
    

    def solve(self):
        n = self.size
        o = NthRootOfUnity(n)
        AT = self.FFT(A, o)
        BT = self.FFT(B, o)
        C = [AT[i]*BT[i] for i in range(n)]
        D = [round((a/n).real) for a in self.FFT(C, o ** -1)]
        while len(D) > 0 and D[-1] == 0:
            del D[-1]
        return D



if __name__ == "__main__":
    A = list(map(int,input().split()))
    B = list(map(int,input().split()))

    C = FPM(A,B)
    ANS = C.solve()
    print(ANS)
