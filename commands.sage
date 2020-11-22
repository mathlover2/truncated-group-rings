load('truncated_group_rings.sage')

from itertools import product as iprod

class TGRExample(object):

    def __init__(self, G, R):

        self.G = G

        self.R = R

        self.RG = self.G.algebra(self.R)

        self.I = TruncatingIdeal(self.RG)

        self.Iq = TruncatingIdeal_quick(self.RG)

        self.RGt = QuotientRing(self.RG, self.I)

        self.RGtq = QuotientRing(self.RG, self.Iq)

        self.fs = self.RG.gens()

        self.fts = self.RGt.gens()

        self.sigma = self.I.gen()

        

Qw = CyclotomicField(3)

Qi = CyclotomicField(4)

QC3 = TGRExample(AbelianGroup([3]), QQ)

QC3C2 = TGRExample(AbelianGroup([3,2]), QQ)

QC3C3 = TGRExample(AbelianGroup([3,3]), QQ)

QwC3C3 = TGRExample(AbelianGroup([3,3]), Qw)

QwC3C2 = TGRExample(AbelianGroup([3,2]), Qw)

QiC4C4 = TGRExample(AbelianGroup([4,4]), Qi)





w = Qw.gen()

f, g = QwC3C3.fts

basis = [(1+f+f^2) * (1/3), (1+g+g^2) * (1/3), (1+f^2*g+f*g^2) * (1/3), (1+f*g+f^2*g^2) * (1/3)]

for units in iprod([1,g,g^2,-1,-g,-g^2], [1,f,f^2,-1,-f,-f^2],[1,f,f^2,-1,-f,-f^2],[1,f,f^2,-1,-f,-f^2] ):
    z = (sum(map(lambda x,y: x*y, units, basis)))

    if all(c in ZZ for c in QwC3C3.RGt.lift(z).coefficients()):

        show(z)



f, g = QiC4C4.fts

ii = Qi.gen()

z1 = (1 + ii*f - f^2 - ii*f^3)*(1 + ii*g - g^2 - ii*g^3)+ (1 - ii*f - f^2 + ii*f^3)*(1 - ii*g - g^2 + ii*g^3)

z2 = (1 - ii*f - f^2 + ii*f^3)*(1 + ii*g - g^2 - ii*g^3) + (1 + ii*f - f^2 - ii*f^3)*(1 - ii*g - g^2 + ii*g^3)

basis = [(1+f+f^2+f^3)*(1-g+g^2-g^3) * (1/16),

         (1+f+f^2+f^3)*(1-g^2) * (1/8),

         (1+g+g^2+g^3)*(1-f^2) * (1/8),

         (1-g+g^2-g^3)*(1-f^2) * (1/8),

         z1 * (1/16),

         z2 * (1/16),

         (1-f+f^2-f^3)*(1+g+g^2+g^3) * (1/16),

         (1-f+f^2-f^3)*(1-g+g^2-g^3) * (1/16),

         (1-f+f^2-f^3)*(1-g^2) * (1/8)]

for units in iprod([1,-1], [1,g,-1,-g], [1,f,-1,-f], [1,f,-1,-f],

                   [1,g,-1,-g], [1,g,-1,-g], [1,-1],

                   [1,-1], [1,g,-1,-g]):

    z = (sum(map(lambda x,y: x*y, units, basis)))

    if all(c.is_integer() for c in QiC4C4.RGt.lift(z).coefficients()):

        show(z)
