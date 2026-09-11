"""Independent arithmetic/algebra checks for the revision notes.

Run: python docs/check_revision_examples.py
Standard library only. These checks are not theorem proofs, coverage certification,
or tests of the research pipeline. The original research code is unchanged.
"""
from __future__ import annotations
import math
import unittest
from collections.abc import Sequence


def transition(theta: float, mu: float, sigma2: float, x: float,
               gap: float) -> tuple[float, float]:
    """Exact OU conditional moments at a fixed or exogenous positive gap."""
    if not all(math.isfinite(z) for z in (theta, mu, sigma2, x, gap)):
        raise ValueError('All inputs must be finite')
    if min(theta, sigma2, gap) <= 0:
        raise ValueError('theta, sigma2 and gap must be positive')
    return (mu + math.exp(-theta*gap)*(x-mu),
            sigma2*(-math.expm1(-2*theta*gap))/(2*theta))


def nll(theta: float, mu: float, sigma2: float,
        x: Sequence[float], gaps: Sequence[float]) -> float:
    """Direct conditional likelihood, separate from the profiling algebra."""
    total = 0.0
    for i, gap in enumerate(gaps):
        m, v = transition(theta, mu, sigma2, x[i], gap)
        total += .5*(math.log(2*math.pi*v)+(x[i+1]-m)**2/v)
    return total


def profile(theta: float, x: Sequence[float], gaps: Sequence[float]
            ) -> tuple[float, float, float]:
    """Profile free mu and sigma squared at fixed theta; no optimization here."""
    if not math.isfinite(theta) or theta <= 0:
        raise ValueError('theta must be finite and positive')
    if not gaps or len(x) != len(gaps)+1:
        raise ValueError('Need n gaps and n+1 observations')
    if not all(math.isfinite(y) for y in x):
        raise ValueError('Observations must be finite')
    if not all(math.isfinite(g) and g > 0 for g in gaps):
        raise ValueError('Gaps must be finite and positive')
    b = [-math.expm1(-theta*g) for g in gaps]
    c = [-math.expm1(-2*theta*g)/(2*theta) for g in gaps]
    d = [x[i+1]-math.exp(-theta*g)*x[i] for i,g in enumerate(gaps)]
    denominator = sum(bi*bi/ci for bi,ci in zip(b,c))
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError('Degenerate profile weights')
    mu = sum(bi*di/ci for bi,di,ci in zip(b,d,c))/denominator
    sigma2 = sum((di-bi*mu)**2/ci for bi,di,ci in zip(b,d,c))/len(gaps)
    if not math.isfinite(sigma2) or sigma2 <= 0:
        raise ValueError('Degenerate or non-finite profiled variance')
    objective = .5*(len(gaps)*(math.log(2*math.pi)+1+math.log(sigma2))
                   + sum(math.log(ci) for ci in c))
    return mu, sigma2, objective


class RevisionExamples(unittest.TestCase):
    def test_coin_likelihood(self):
        likelihood = math.comb(10,7)*.7**7*.3**3
        self.assertAlmostEqual(likelihood,.266827932)
        self.assertEqual(round(math.log(likelihood),4),-1.3212)
        self.assertAlmostEqual(7/.7-3/.3,0)

    def test_walk_moments(self):
        m = 2/4-3/4
        v = 4/4+3/4-m*m
        self.assertEqual((m,v,100*m,100*v),(-.25,27/16,-25,168.75))

    def test_walk_overlap(self):
        self.assertFalse(set(range(3,6)) & set(range(6,9)))
        self.assertEqual(len(set(range(3,6)) & set(range(5,7))),1)

    def test_ruin_recursion_and_boundaries(self):
        for p in (.3,.5,.8):
            q=1-p
            f=lambda m: m/10 if p==.5 else (1-(q/p)**m)/(1-(q/p)**10)
            self.assertAlmostEqual(f(0),0)
            self.assertAlmostEqual(f(10),1)
            for m in range(1,10):
                self.assertAlmostEqual(f(m),p*f(m+1)+q*f(m-1))

    def test_stationary_chain(self):
        pi=(4/7,3/7)
        self.assertAlmostEqual(.7*pi[0]+.4*pi[1],pi[0])
        self.assertAlmostEqual(.3*pi[0]+.6*pi[1],pi[1])

    def test_endpoint_vs_route(self):
        self.assertAlmostEqual(.7**2+.3*.4,.61)
        self.assertAlmostEqual(.7**2,.49)

    def test_weather_bridge(self):
        self.assertAlmostEqual(.2*.6/(.8*.2+.2*.6),3/7)

    def test_ar_moments(self):
        self.assertEqual(25+.5*(29-25),27)
        self.assertEqual(25+.5**2*(29-25),26)
        self.assertEqual(4*(1-.5**4)/(1-.5**2),5)
        self.assertEqual(4/(1-.5**2),16/3)

    def test_ar_regression(self):
        u,y=(0,2,1),(2,1,3)
        um,ym=sum(u)/3,sum(y)/3
        b=sum((a-um)*(b-ym) for a,b in zip(u,y))/sum((a-um)**2 for a in u)
        a=ym-b*um
        rss=sum((yi-a-b*ui)**2 for ui,yi in zip(u,y))
        self.assertEqual((b,a,rss,rss/3),(-.5,2.5,1.5,.5))
        self.assertAlmostEqual(a/(1-b),5/3)

    def test_path_events(self):
        paths=[(a,a+b) for a in (-1,1) for b in (-1,1)]
        endpoint={p for p in paths if p[1]==0}
        positive={p for p in paths if min(p)>=0}
        self.assertEqual((len(endpoint),len(positive),len(endpoint & positive)),(2,2,1))

    def test_interpolated_variance(self):
        self.assertEqual((0+.5**2)/4,1/16)
        self.assertNotEqual((0+.5**2)/4,1/8)

    def test_brownian_moments(self):
        self.assertEqual(min(2,5),2)
        self.assertAlmostEqual(2/math.sqrt(2*5),math.sqrt(2/5))
        self.assertEqual(5-2,3)

    def test_brownian_overlap(self):
        self.assertEqual(min(3,4)-min(3,2)-min(1,4)+min(1,2),1)

    def test_constant_sde_and_weighted_noise(self):
        t=2
        self.assertEqual((1+3*t,2**2*t),(7,8))
        self.assertEqual(2**2*1+(-1)**2*2,6)

    def test_ou_table(self):
        expected={.1:(28.6193,.3625),.5:(27.4261,1.2642),1:(26.4715,1.7293),3:(25.1991,1.9950)}
        for gap,pair in expected.items():
            self.assertEqual(tuple(round(v,4) for v in transition(1,25,4,29,gap)),pair)

    def test_semigroup_and_stationarity(self):
        for theta in (.02,1,5):
            for h,k in ((.001,.4),(.5,1.2),(3,8)):
                m1,v1=transition(theta,25,4,29,h)
                m2,v2=transition(theta,25,4,m1,k)
                m,v=transition(theta,25,4,29,h+k)
                self.assertAlmostEqual(m2,m)
                self.assertAlmostEqual(math.exp(-2*theta*k)*v1+v2,v)
                q=4/(2*theta)
                self.assertAlmostEqual(math.exp(-2*theta*h)*q+v1,q)

    def test_ar_to_ou_and_short_gap(self):
        theta,sigma2=math.log(2),(32/3)*math.log(2)
        m,v=transition(theta,25,sigma2,29,1)
        self.assertAlmostEqual(m,27)
        self.assertAlmostEqual(v,4)
        _,v=transition(1,25,4,29,1e-10)
        self.assertAlmostEqual(v/(4e-10),1,places=9)

    def test_profile_objective_and_minimum(self):
        x=(2,2.3,1.8,2.5,2,1.7,2.1)
        gaps=(.1,.7,.25,1.2,.4,.9)
        for theta in (.05,.8,3):
            mu,sigma2,value=profile(theta,x,gaps)
            self.assertAlmostEqual(value,nll(theta,mu,sigma2,x,gaps),places=11)
            for delta in (-.1,.1):
                self.assertGreater(nll(theta,mu+delta,sigma2,x,gaps),value)
            for factor in (.8,1.2):
                self.assertGreater(nll(theta,mu,sigma2*factor,x,gaps),value)

    def test_profile_invalid_and_degenerate(self):
        for theta,x,g in ((0,[1,2],[1]),(1,[1,2],[0]),(1,[1,math.nan],[1]),(1,[1,2],[]),(1,[0,0,0],[1,1])):
            with self.assertRaises(ValueError):
                profile(theta,x,g)

    def test_exponential_gap_target(self):
        self.assertAlmostEqual(math.log(1+.5)/.5,.8109302162163288)
        self.assertLess(math.log(1+.5)/.5,1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
