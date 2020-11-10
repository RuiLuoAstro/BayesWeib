from utils import *
import numpy as np
import scipy.special as sp
import pymultinest
import warnings
import sys
import argparse

dis = StatDis()
ld = LoadDat()

d2hr = 24.

def lnlik(vpar):
    try:
        p1 = dis.ProbNN(ts1_o1, tne_o1, delta_o1, vpar[0], vpar[1])
        p2 = dis.ProbN0(tobs_o2, vpar[0], vpar[1])
        p3 = dis.ProbNN(ts1_o3, tne_o3, delta_o3, vpar[0], vpar[1])
        p4 = dis.ProbNN(ts1_o4, tne_o4, delta_o4, vpar[0], vpar[1])
        loglik = np.log(p1) + np.log(p2) + np.log(p3) + np.log(p4)
        return loglik
    except:
        print 'Numerical error: @', vpar
        return -1e99

def prior(cube, ndim, nparams):
    for i in range(ndim):
        cube[i] = vpar_range[i,0]+cube[i] *(vpar_range[i,1]-vpar_range[i,0])
        cube[i] = np.power(10, cube[i])

def loglikfunc(cube, ndim, nparams):
    '''
    The log likelihood function. This function has to be called "LogLikelihood".

    Args:
        cube (:class:`numpy.ndarray`): an array of parameter values.

    Returns:
        float: the log likelihood value.
    '''
    cube2 = np.zeros(ndim)
    for i in range(0,ndim):
        cube2[i] = cube[i]
    res = lnlik(cube2)
    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Weibull distribution analysis')
    parser.add_argument('-tim', action='store', dest='ftim', type=str, help='Input the file of Arrival Times')
    parser.add_argument('-obs', action='store', dest='fobs', type=str, help='Input the file of Observations')
    parser.add_argument('-out', action='store', dest='fout', type=str, help='Save the output file')

    args = parser.parse_args()
    ftim = args.ftim
    fobs = args.fobs
    fout = args.fout

    rptinfo = ld.LoadTim(ftim)
    obsinfo = ld.LoadObs(fobs)
    
    toa = rptinfo['ToA']
    obs = rptinfo['Obs']
    toa_o1 = toa[obs=='1st']
    toa_o3 = toa[obs=='3rd']
    toa_o4 = toa[obs=='4th']
    tobs = obsinfo['T_end'] - obsinfo['T_start']
    tobs_o2 = tobs[1]
    delta_o1 = toa_o1[1:] - toa_o1[:-1]
    delta_o3 = toa_o3[1:] - toa_o3[:-1]
    delta_o4 = toa_o4[1:] - toa_o4[:-1]
    ts1_o1 = toa_o1[0] - obsinfo['T_start'][0]
    tne_o1 = obsinfo['T_end'][0] - toa_o1[-1]
    ts1_o3 = toa_o3[0] - obsinfo['T_start'][2]
    tne_o3 = obsinfo['T_end'][2] - toa_o3[-1]
    ts1_o4 = toa_o4[0] - obsinfo['T_start'][3]
    tne_o4 = obsinfo['T_end'][3] - toa_o4[-1]

    # Transform the unit from day to hour
    tobs_o2 = tobs_o2*d2hr
    delta_o1 = delta_o1*d2hr
    delta_o3 = delta_o3*d2hr
    delta_o4 = delta_o4*d2hr
    ts1_o1 = ts1_o1*d2hr
    tne_o1 = tne_o1*d2hr
    ts1_o3 = ts1_o3*d2hr
    tne_o3 = tne_o3*d2hr
    ts1_o4 = ts1_o4*d2hr
    tne_o4 = tne_o4*d2hr

    vpara=np.array([-1, -2])
    vparb=np.array([1, 1])
    vpar_range=np.dstack((vpara.transpose(),vparb.transpose()))[0,:,:]

    print '-----------------par range------------------'
    print vpar_range
    print 
    print '-------------marginal likelihood------------'
    print loglikfunc(vpara, len(vpara),len(vpara))
    print 
    print "The Nest sampler is running ..."
    print '-------------sampling process---------------'
    # run MultiNest
    pymultinest.run(loglikfunc, prior, len(vpara),
                    importance_nested_sampling = False,
                    resume = False,
                    verbose = True,
                    sampling_efficiency = 'model',
                    n_live_points = 5000,
                    outputfiles_basename='nest_out/'+fout)
