import numpy as np
import scipy.special as sp
from scipy.interpolate import interp1d

class MathFunc():
    def gammainc(self, alpha, x):
        '''
        Incomplete Gamma function including all of cases.
        Args:
        alpha: Power-law index argument
        x: Lower limit of integral, x must >= 0
        '''
        if alpha==0:
            return -sp.expi(-x)
        elif alpha<0:
            return (gammainc(alpha+1,x)-np.power(x, alpha)*np.exp(-x))/alpha
        else:
            return sp.gammaincc(alpha,x)*sp.gamma(alpha)

    def Sampling1D(self, x, y, x1, x2, n):
        '''
        Monte Carlo sampling for 1-D parameter
        x: Independent variable range
        y: Function range
        x1: Lower border of samples in x-space
        x2: Upper border of samples x-space
        n: Sampling number
        '''
        nt = 0
        res = np.array([])
        while (nt<n):
            fuc = interp1d(x, y / np.max(y))
            vx = np.random.uniform(x1, x2, n-nt)
            vy = np.random.uniform(0, 1, n-nt)
            res = np.append(res, vx[vy <= fuc(vx)])
            nt = len(res)
        return res

    def SamplingND(self, fuc, par_range, maxv_ori, n):
        '''
        Monte Carlo sampling for N-D parameter
        func: Function range, array
        par_range: Independent variable ranges, array
        maxv_ori: Border ranges of samples, array
        n: Sampling number
        '''
        nt = 0
        maxv = maxv_ori
        res = np.array([])
        npar, m = par_range.shape
        res = res.reshape((0, npar))
        while (nt < n):
            vpar = np.random.uniform(0, 1, (n-nt, npar))
            for i in range(npar):
                lv = par_range[i,0]
                rv = par_range[i,1]
                vpar[:,i] = vpar[:,i] * (rv-lv) + lv
            
            vy = np.random.uniform(0, 1, n-nt)*maxv
            fv = fuc(vpar)
            if (np.max(fv) > maxv):
                maxv = np.max(fv)*1.5
                nt = 0
                res = np.array([])
                res = res.reshape((0, npar))
            else:
                res = np.vstack( (res, vpar[vy <= fv,:]))
                nt, m = res.shape
            print nt
        return res

class StatDis():
    def __init__(self):
        self.mf = MathFunc()

    def Weibull(self, delta, r, k):
        '''
        The Probability density function of Weibull distribution.
        Args:
            delta: Time interval of adjacent events in untis of sec.
            r: Event rate in units of hr^{-1},
            k: The shape parameter of Weibull distribution, for example, for k=1, Weibull reduces to exponential distribution that follow the Poisson process.
        '''
        lamda = delta*r*sp.gamma(1+1./k)
        pd = k/delta*np.power(lamda, k)*np.exp(-np.power(lamda, k))
        return pd

    def logWeibull(self, delta, r, k):
        '''
        Logarithmic probablity density of Weibull distribution.
        Args:
            delta: Time interval of adjacent events in untis of sec.
            r: Event rate in units of hr^{-1},
            k: The shape parameter of Weibull distribution, for example, for k=1, the Weibull distribution reduces to exponential one that follow the Poissonian process.
        '''
        par_a = np.log(k)+(k-1)*np.log(delta)+k*np.log(r)+sp.gammaln(1+1./k)
        par_b = np.power(delta*r*sp.gamma(1+1./k), k)
        logpd = par_a - par_b
        return logpd

    def CDF(self, delta, r, k):
        '''
        Cumulative distribution function of Weibull distribution.
        '''
        lamda = delta*r*sp.gamma(1+1./k)
        return np.exp(-np.power(lamda, k))

    def ProbNN(self, ts1, tne, delta, r, k):
        '''
        Probability of N events detected in tobs observing time.
        Args:
            ts1: Time interval from the start time to the first event
            tne: Time interval from the last event to the end time
        '''
        w = self.Weibull(delta, r, k)
        res = r*self.CDF(ts1, r, k)*self.CDF(tne, r, k)*np.prod(w)
        return res

    def ProbN1(self, ts1, t1e, r, k):
        '''
        Probability of one event detected in tobs observing time.
        Args:
            ts1: Time interval from the start time to this event
            t1e: Time interval from this event to the end time
        '''
        res = r*self.CDF(ts1, r, k)*self.CDF(t1e, r, k)
        return res

    def ProbN0(self, tobs, r, k):
        '''
        Probability of zero event detected in tobs observing time.
        Args:
            tobs: A single observation duration
        '''
        lamda = tobs*r*sp.gamma(1+1./k)
        res = self.mf.gammainc(1./k, np.power(lamda, k))/k/sp.gamma(1+1./k)
        return res


class LoadDat():
    def LoadTim(self, fin):
        cat = np.loadtxt(fin, dtype='string')
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0,i]] = cat[1:, i]
        cat2['Nb'] = np.array(cat2['Nb'], dtype=int)
        cat2['ToA'] = np.array(cat2['ToA'], dtype=float)
        cat2['Obs'] = np.array(cat2['Obs'])
        return cat2

    def LoadObs(self, fin):
        cat = np.loadtxt(fin, dtype='string')
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0,i]] = cat[1:, i]
        cat2['Obs'] = np.array(cat2['Obs'])
        cat2['T_start'] = np.array(cat2['T_start'], dtype=float)
        cat2['T_end'] = np.array(cat2['T_end'], dtype=float)
        return cat2