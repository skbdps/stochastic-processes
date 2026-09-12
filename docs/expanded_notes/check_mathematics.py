"""Independent numerical/algebraic checks for the teaching notes.

This editor-side suite uses NumPy and SciPy. It checks worked examples and
identities; it does not certify general theorems, research results or all
possible numerical inputs. The accompanying example module is standard-library only.
"""
from __future__ import annotations
import itertools
import json
import math
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from ou_likelihood_examples import conditional_nll, profile_at_theta, transition_coefficients, validate_data


def ar_fit(x):
    x=np.asarray(x,dtype=float);u=x[:-1];y=x[1:]
    b=np.sum((u-u.mean())*(y-y.mean()))/np.sum((u-u.mean())**2)
    a=y.mean()-b*u.mean()
    v=np.mean((y-a-b*u)**2)
    return a,b,v


def ou_moments(theta,mu,sigma,x,h):
    phi,b,c=transition_coefficients(theta,h)
    return phi*x+b*mu,sigma*sigma*c


class TeachingExamples(unittest.TestCase):
    def test_coin_record_and_count(self):
        record=(1,1,0,1,0,1,1,1,0,1)
        self.assertEqual(sum(record),7);self.assertEqual(math.comb(10,7),120)
        self.assertAlmostEqual(math.comb(10,7)*.7**7*.3**3,.266827932)

    def test_coin_candidate_table(self):
        expected={.1:.000008748,.3:.009001692,.5:.1171875,.7:.266827932,.9:.057395628}
        for p,l in expected.items():self.assertAlmostEqual(120*p**7*(1-p)**3,l,places=10)

    def test_coin_score_and_maximum(self):
        self.assertAlmostEqual(7/.7-3/.3,0)
        self.assertGreater(7/.4-3/.6,0);self.assertLess(7/.9-3/.1,0)
        self.assertEqual(120*0**7*1**3,0);self.assertEqual(120*1**7*0**3,0)

    def test_coin_estimator_moments_by_enumeration(self):
        p=.3;n=10
        probs=[math.comb(n,k)*p**k*(1-p)**(n-k) for k in range(n+1)]
        mean=sum((k/n)*v for k,v in enumerate(probs))
        var=sum((k/n-mean)**2*v for k,v in enumerate(probs))
        self.assertAlmostEqual(mean,p);self.assertAlmostEqual(var,p*(1-p)/n)

    def test_gaussian_mle_moments(self):
        x=np.array([1.,2.,5.]);self.assertAlmostEqual(x.mean(),8/3)
        self.assertAlmostEqual(np.mean((x-x.mean())**2),26/9)

    def test_walk_enumeration(self):
        walks=np.cumsum(np.array(list(itertools.product([-1,1],repeat=3))),axis=1)
        self.assertEqual(len(walks),8)
        self.assertAlmostEqual(walks[:,-1].mean(),0);self.assertAlmostEqual(walks[:,-1].var(),3)
        self.assertAlmostEqual(np.mean(walks[:,0]*walks[:,2]),1)

    def test_walk_unequal_steps(self):
        mean=.25*2+.75*(-1);v=.25*4+.75-mean**2
        self.assertEqual(mean,-.25);self.assertEqual(v,27/16)
        self.assertEqual(100*v,168.75);self.assertEqual(40*v,67.5)

    def test_walk_overlap(self):
        a=np.array([0,0,1,1,1,0,0,0]);b=np.array([0,0,0,0,1,1,0,0])
        c=np.array([0,0,0,0,0,1,1,1]);self.assertEqual(a@b,1);self.assertEqual(a@c,0)

    def test_gambler_fair_recursion(self):
        n=100
        for m in range(1,n):self.assertAlmostEqual(m/n,.5*((m+1)/n+(m-1)/n))
        self.assertEqual(30/100,.3)

    def test_gambler_duration(self):
        n=100;u=lambda m:m*(n-m)
        for m in range(1,n):self.assertEqual(u(m),1+.5*u(m+1)+.5*u(m-1))
        self.assertEqual(u(30),2100)

    def test_gambler_biased_recursion(self):
        for p,n in [(2/3,3),(.4,7),(.8,12)]:
            r=(1-p)/p;f=lambda m:(1-r**m)/(1-r**n)
            for m in range(1,n):self.assertAlmostEqual(f(m),p*f(m+1)+(1-p)*f(m-1))
        self.assertAlmostEqual((1-.5)/(1-.5**3),4/7)

    def test_shifted_boundaries(self):
        A,B,m=2,3,0
        self.assertEqual((m+A)/(A+B),2/5);self.assertEqual((m+A)*(B-m),6)

    def test_markov_matrix(self):
        a=np.array([[.7,.4],[.3,.6]])
        np.testing.assert_allclose(a@a,[[.61,.52],[.39,.48]])
        pi=np.array([4/7,3/7]);np.testing.assert_allclose(a@pi,pi)
        np.testing.assert_allclose(a@np.array([.5,.5]),[.55,.45])

    def test_markov_convergence_formula(self):
        a=np.array([[.7,.4],[.3,.6]])
        for n in (0,1,2,7,20):
            val=(np.linalg.matrix_power(a,n)@np.array([1.,0.]))[0]
            self.assertAlmostEqual(val,4/7+.3**n*(1-4/7))

    def test_periodic_and_one_step_examples(self):
        flip=np.array([[0.,1.],[1.,0.]])
        np.testing.assert_allclose(flip@flip,np.eye(2))
        pi=np.array([.2,.8]);a=np.column_stack([pi,pi]);np.testing.assert_allclose(a@np.array([1.,0.]),pi)

    def test_record_vs_endpoint(self):
        self.assertAlmostEqual(.5*.7*.3*.4,.042)
        self.assertAlmostEqual(.7*.3*.4,.084)
        self.assertAlmostEqual(.7**2+.3*.4,.61)

    def test_transition_counts(self):
        x=[1,1,2,1,2,2,1,1,1,2,1]
        counts={(i,j):sum(a==i and b==j for a,b in zip(x[:-1],x[1:])) for i in (1,2) for j in (1,2)}
        self.assertEqual(counts,{(1,1):3,(1,2):3,(2,1):3,(2,2):1})
        self.assertEqual(counts[(1,1)]/(counts[(1,1)]+counts[(1,2)]),.5)

    def test_weather_bridge(self):
        self.assertAlmostEqual(.2*.6/(.8*.2+.2*.6),3/7)

    def test_ar_gaussian_record_density(self):
        density=(1/math.sqrt(8*math.pi))**2*math.exp(-(.5**2)/(2*4))
        self.assertAlmostEqual(density,math.exp(-1/32)/(8*math.pi))

    def test_ar_moment_geometric_sum(self):
        for phi in (.5,-.5,.9):
            for k in (1,2,5):
                self.assertAlmostEqual(sum(4*phi**(2*r) for r in range(k)),4*(1-phi**(2*k))/(1-phi**2))
        self.assertEqual(25+.5**2*4,26)
        self.assertEqual(4*(1+.5**2),5)

    def test_ar_stationary_variance(self):
        for phi in (.5,-.5,.9):
            q=4/(1-phi**2);self.assertAlmostEqual(phi**2*q+4,q)
        self.assertEqual(4/(1-.25),16/3)

    def test_ar_positive_fit(self):
        a,b,v=ar_fit([0,1,1,2])
        self.assertAlmostEqual(a,1);self.assertAlmostEqual(b,.5);self.assertAlmostEqual(v,1/6)
        self.assertAlmostEqual(a/(1-b),2);self.assertAlmostEqual(v/(1-b*b),2/9)

    def test_ar_negative_fit(self):
        a,b,v=ar_fit([0,2,1,3])
        self.assertAlmostEqual(a,2.5);self.assertAlmostEqual(b,-.5);self.assertAlmostEqual(v,.5)
        self.assertAlmostEqual(a/(1-b),5/3)

    def test_ar_constrained_refit(self):
        x=np.array([0.,2.,1.,3.]);u=x[:-1];y=x[1:]
        ss=lambda b:sum((y-y.mean()-b*(u-u.mean()))**2)
        self.assertTrue(all(ss(0)<=ss(b)+1e-12 for b in np.linspace(0,.9,50)))

    def test_clock_moment_calibration(self):
        delta=.01;n=100;a=2;sigma=3
        self.assertAlmostEqual(n*a*delta,2)
        self.assertAlmostEqual(n*(sigma*math.sqrt(delta))**2,9)
        self.assertAlmostEqual(n*(sigma*delta)**2,.09)

    def test_interpolation_variance(self):
        n=4;r=.5;k=0
        self.assertEqual((k+r*r)/n,1/16);self.assertEqual((k+r)/n,1/8)
        self.assertEqual(((k+r)-(k+r*r))/n,1/(4*n))

    def test_brownian_covariance_matrix(self):
        t=np.array([1.,3.,4.]);c=np.minimum.outer(t,t)
        np.testing.assert_allclose(c,[[1,1,1],[1,3,3],[1,3,4]])
        self.assertTrue(np.all(np.linalg.eigvalsh(c)>0))

    def test_brownian_overlap(self):
        cov=min(3,4)-min(3,2)-min(1,4)+min(1,2)
        self.assertEqual(cov,1);self.assertEqual(cov/math.sqrt(2*2),.5)

    def test_gaussian_fourth_moment(self):
        ans=quad(lambda z:z**4*math.exp(-z*z/2)/math.sqrt(2*math.pi),-np.inf,np.inf)[0]
        self.assertAlmostEqual(ans,3,places=8)

    def test_quadratic_variation_variance(self):
        h=np.array([.1,.3,.2,.4]);var=2*np.sum(h*h)
        self.assertLessEqual(var,2*h.sum()*h.max()+1e-12)

    def test_wiener_weights_and_covariance(self):
        f=np.array([2.,-1.,-1.]);g=np.array([1.,1.,0.])
        self.assertEqual(f@f,6);self.assertEqual(g@g,2);self.assertEqual(f@g,1)
        self.assertEqual((f+g)@(f+g),10)

    def test_wiener_smooth_integrals(self):
        for t in (.2,1.,3.):
            self.assertAlmostEqual(quad(lambda s:s*s,0,t)[0],t**3/3)
            for theta in (.1,1.,3.):
                value=quad(lambda s:math.exp(-2*theta*(t-s)),0,t)[0]
                self.assertAlmostEqual(value,-math.expm1(-2*theta*t)/(2*theta))

    def test_integration_by_parts_identity(self):
        # A smooth deterministic path checks the finite-partition algebra.
        theta=.7;t=2.;f=lambda s:math.exp(-theta*(t-s))
        lhs=quad(lambda s:f(s)*2*s,0,t)[0]
        rhs=t*t-quad(lambda s:theta*f(s)*s*s,0,t)[0]
        self.assertAlmostEqual(lhs,rhs)

    def test_ou_pathwise_integral_equation(self):
        # The noise-subtracted construction solves the integral equation for
        # every continuous driving path, so use sin(s) as a quadrature check.
        theta=1.4;mu=2.;sigma=.7;x0=3.
        def x(t):
            return mu+math.exp(-theta*t)*(x0-mu)+sigma*(math.sin(t)-theta*quad(lambda s:math.exp(-theta*(t-s))*math.sin(s),0,t)[0])
        for t in (.2,1.,2.):
            rhs=x0+theta*quad(lambda s:mu-x(s),0,t)[0]+sigma*math.sin(t)
            self.assertAlmostEqual(x(t),rhs,places=9)

    def test_ou_transition_table(self):
        expected={.1:(28.6193,.3625),.5:(27.4261,1.2642),1.:(26.4715,1.7293),3.:(25.1991,1.9950)}
        for dt,(m,v) in expected.items():
            got=ou_moments(1,25,2,29,dt)
            self.assertEqual(round(got[0],4),m);self.assertEqual(round(got[1],4),v)

    def test_ou_semigroup(self):
        for theta in (.02,1.,7.):
            for h,k in ((.001,.3),(.5,1.2),(3.,8.)):
                m1,v1=ou_moments(theta,25,2,29,h)
                m2,v2=ou_moments(theta,25,2,m1,k)
                m,v=ou_moments(theta,25,2,29,h+k)
                self.assertAlmostEqual(m,m2,places=10)
                self.assertAlmostEqual(v,math.exp(-2*theta*k)*v1+v2,places=10)

    def test_ou_stationarity_and_small_gap(self):
        for theta in (.1,1.,3.):
            q=4/(2*theta)
            for h in (.001,.5,10.):
                phi,b,c=transition_coefficients(theta,h)
                self.assertAlmostEqual(phi*phi*q+4*c,q)
        _,v=ou_moments(1,25,2,29,1e-10)
        self.assertAlmostEqual(v/(4e-10),1,places=8)

    def test_ou_fixed_start_covariance(self):
        theta=.7;q=2.
        for s,t in ((.2,.8),(1.,3.),(0.,2.)):
            a=math.exp(-theta*(t-s))*q*(1-math.exp(-2*theta*s))
            b=q*(math.exp(-theta*(t-s))-math.exp(-theta*(t+s)))
            self.assertAlmostEqual(a,b)

    def test_ou_parameter_conversion(self):
        theta=math.log(2);sigma2=(32/3)*math.log(2)
        m,v=ou_moments(theta,25,math.sqrt(sigma2),29,1)
        self.assertAlmostEqual(m,27);self.assertAlmostEqual(v,4)
        self.assertAlmostEqual(math.exp(-theta*(math.log(2)/theta)),.5)

    def test_ch13_candidate_values(self):
        x=(29,27.5,26,25.8);dt=(.2,.8,.3)
        expected={.5:3.978404677316483,1:3.1931010128810344,2:2.461343441684951}
        for th,value in expected.items():self.assertAlmostEqual(conditional_nll(th,25,2,x,dt),value)

    def test_profile_table_and_direct_objective(self):
        x=(29,27.5,26,25.8);dt=(.2,.8,.3)
        exp={.5:(21.819131632369473,1.498570868792321,3.0362759028433417),1:(24.240779169404064,.9906802821682064,2.137539544990499),2:(25.4040127609971,.2900227368374021,-.17762806104800655)}
        for th,(mu,sig2,nll) in exp.items():
            fit=profile_at_theta(th,x,dt)
            self.assertAlmostEqual(fit.mu,mu);self.assertAlmostEqual(fit.sigma2,sig2);self.assertAlmostEqual(fit.nll,nll)
            self.assertAlmostEqual(conditional_nll(th,fit.mu,math.sqrt(fit.sigma2),x,dt),fit.nll,places=9)

    def test_profile_against_independent_optimizer(self):
        x=(2,2.3,1.8,2.5,2,1.7,2.1);dt=(.1,.7,.25,1.2,.4,.9)
        for theta in (.05,.8,3.):
            fit=profile_at_theta(theta,x,dt)
            objective=lambda z:conditional_nll(theta,float(z[0]),math.exp(float(z[1])),x,dt)
            result=minimize(objective,[fit.mu+.2,math.log(math.sqrt(fit.sigma2))+.1],method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-10})
            self.assertTrue(result.success)
            self.assertAlmostEqual(result.fun,fit.nll,places=8)
            self.assertAlmostEqual(result.x[0],fit.mu,places=5)

    def test_profile_local_minima_of_nuisance_parameters(self):
        x=(2,2.3,1.8,2.5,2);dt=(.1,.7,.25,1.2)
        fit=profile_at_theta(.8,x,dt)
        for delta in (-.1,.1):
            self.assertGreater(conditional_nll(.8,fit.mu+delta,math.sqrt(fit.sigma2),x,dt),fit.nll)
        for factor in (.8,1.2):
            self.assertGreater(conditional_nll(.8,fit.mu,math.sqrt(fit.sigma2*factor),x,dt),fit.nll)

    def test_exponential_gap_target(self):
        theta=1.;d=.5
        val=quad(lambda u:math.exp(-theta*u)*math.exp(-u/d)/d,0,np.inf)[0]
        self.assertAlmostEqual(val,1/(1+theta*d))
        target=-math.log(val)/d
        self.assertAlmostEqual(target,.8109302162163288)
        self.assertLess(target,theta)

    def test_jensen_two_point_gaps(self):
        theta=1.;a=.1;b=.9;mean=(a+b)/2
        laplace=.5*(math.exp(-theta*a)+math.exp(-theta*b))
        self.assertGreater(laplace,math.exp(-theta*mean))
        self.assertLess(-math.log(laplace)/mean,theta)

    def test_invalid_input_rejection(self):
        for x,g in (([1,2],[0]),([1,2],[-1]),([1,2],[]),([1,math.nan],[1]),([1,2],[math.inf])):
            with self.assertRaises(ValueError):validate_data(x,g)
        for theta in (0,-1,math.inf,math.nan):
            with self.assertRaises(ValueError):transition_coefficients(theta,.2)
        with self.assertRaises(ValueError):profile_at_theta(1,[0,0,0],[1,1])

    def test_sde_example_moments(self):
        t=2.;self.assertEqual(1+3*t,7);self.assertEqual(4*t,8)
        self.assertEqual(4+t*t+t,10);self.assertEqual(9*t,18)


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(TeachingExamples)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    report={'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'successful':result.wasSuccessful(),'scope':'Worked-example arithmetic and selected algebraic/numerical identities; not general theorem or peer-review certification.'}
    (Path(__file__).parent/'validation'/'mathematics_report.json').write_text(json.dumps(report,indent=2))
    raise SystemExit(0 if result.wasSuccessful() else 1)
