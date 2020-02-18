#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
import numpy as np

def basicBootCImean(data, nresamp=1000, alpha=0.95):
    '''
    Return (basic) bootstrapped confidenced interval for mean on data
    
    data contains a list of numeric values.  This function uses basic,
    non-parametric bootstrapping (i.e., no assumptions about the underlying
    distribution of these values) to compute a confidence interval
    (with alpha level of confidence, i.e., the statistic lies within the CI
    with probability alpha) for the mean statistic on these
    values.  nsamp bootstrap resamplings are used.
    '''
        
    thetaHat = float(np.mean(data))
    thetaTildeList = []
    
    if nresamp <= 0:
        raise ValueError('nresamp must be positive!')
        
    for _ in range(nresamp):
        resamp_data = list(np.random.choice(data, len(data), replace=True))
        thetaTilde = float(np.mean(resamp_data))
        thetaTildeList.append(thetaTilde)
    
    thetaTildeList.sort()
    
    indLower = int(np.round(nresamp * ((1.0 - alpha) / 2))) - 1 # 1 - alpha is, e.g., 5%
    indUpper = nresamp - 1 - indLower - 1
    if indLower < 0:
        indLower = 0
        indUpper = nresamp - 1
        raise RuntimeError('Too few resamplings to get desired confidence level!')
    
    lowerCI = 2 * thetaHat - thetaTildeList[indUpper]
    upperCI = 2 * thetaHat - thetaTildeList[indLower]
        
    return lowerCI, upperCI


def studentizedBootCImean(data, nresamp=1000, alpha=0.95, innerResamp=1000):
    '''
    Return (studentized) bootstrapped confidenced interval for mean on data
    
    data contains a list of numeric values.  This function uses studentized,
    non-parametric bootstrapping (i.e., no assumptions about the underlying
    distribution of these values) to compute a confidence interval
    (with alpha level of confidence, i.e., the statistic lies within the CI
    with probability alpha) for the mean statistic on these
    values.  nsamp bootstrap resamplings are used.
    
    Uses a divisor of len(data) to compute the sample std. dev.
    '''
        
    thetaHat = float(np.mean(data))
    thetaTildeList = []
    sTildeList = []
    
    if nresamp <= 0:
        raise ValueError('nresamp must be positive!')
    
    
    for b in range(nresamp):
        resamp_data = list(np.random.choice(data, len(data), replace=True))
        thetaTilde_b = float(np.mean(resamp_data))
        thetaTildeList.append(thetaTilde_b)
        
        thetaTildeInnerList = []
        for _ in range(innerResamp):
            inner_data = list(np.random.choice(resamp_data, len(data), replace=True))
            thetaTildeInner = float(np.mean(inner_data))
            thetaTildeInnerList.append(thetaTildeInner)
            
        sb = float(np.std(thetaTildeInnerList))
        sTildeList.append(sb)
    
    s = float(np.std(thetaTildeList))
    
    fracList = []
    for b in range(nresamp):
        fracList.append((thetaTildeList[b] - thetaHat) / sTildeList[b])
    fracList.sort()
    
    indLower = int(np.round(nresamp * ((1.0 - alpha) / 2))) - 1 # 1 - alpha is, e.g., 5%
    indUpper = nresamp - 1 - indLower - 1
    if indUpper < indLower:
        raise RuntimeError('indUpper is less than indLower!')
    if indLower < 0:
        indLower = 0
        indUpper = nresamp - 1
        raise RuntimeError('Too few resamplings to get desired confidence level!')
            
    lowerCI = thetaHat - s * fracList[indUpper]
    upperCI = thetaHat - s * fracList[indLower]
        
    return lowerCI, upperCI


def bootCImean(data, nresamp=999, alpha=0.95):
    '''
    Return (studentized) bootstrapped confidenced interval for mean on data
    
    data contains a list of numeric values.  This function uses studentized,
    non-parametric bootstrapping (i.e., no assumptions about the underlying
    distribution of these values) to compute a confidence interval
    (with alpha level of confidence) for the mean statistic on these
    values.  nsamp bootstrap resamplings are used.
    '''
    
    if nresamp <= 0:
        raise ValueError('nresamp must be positive!')
    
    meanlist = []
    sigmalist = []
    tlist = []
    
    sampmean = float(np.mean(data))
    templist= []
    for y in data:
        templist.append((y - sampmean) ** 2)
    sampsigma = float(np.sqrt(float(np.sum(templist))/(len(data)* (len(data) - 1))))
    
    tindUpper = int(np.round((nresamp + 1) * ((1.0 - alpha) / 2))) - 1
    tindLower = (nresamp + 1 ) - 1 - tindUpper
    
    for _ in range(nresamp):
        resamp_data = list(np.random.choice(data,len(data),replace=True))
        resamp_mean = np.mean(resamp_data)
        meanlist.append(resamp_mean)
        for y in resamp_data:
            templist.append((y - resamp_mean) ** 2)
        tempsigma = float(np.sqrt(float(np.sum(templist))/(len(resamp_data) * (len(resamp_data) - 1))))
        sigmalist.append(tempsigma)
        tlist.append((resamp_mean - sampmean)/tempsigma)
    
    tlist.sort()
    
    upperCI = sampmean - sampsigma * tlist[tindUpper - 1]
    lowerCI = sampsigma * tlist[tindLower - 1] - sampmean
    
    print('In bootCImean:')
    print('sampsigma =', sampsigma)
    print('selected tlist is:', tlist[0::100] + [tlist[-1]])
    print(sampsigma * tlist[tindUpper - 1], sampsigma * tlist[tindLower - 1])
    
    return lowerCI, upperCI


def basicBootCImedian(data, nresamp=1000, alpha=0.95):
    '''
    Return (basic) bootstrapped confidenced interval for median on data
    
    data contains a list of numeric values.  This function uses basic,
    non-parametric bootstrapping (i.e., no assumptions about the underlying
    distribution of these values) to compute a confidence interval
    (with alpha level of confidence, i.e., the statistic lies within the CI
    with probability alpha) for the median statistic on these
    values.  nsamp bootstrap resamplings are used.
    '''
        
    thetaHat = float(np.median(data))
    thetaTildeList = []
    
    if nresamp <= 0:
        raise ValueError('nresamp must be positive!')
        
    for _ in range(nresamp):
        resamp_data = list(np.random.choice(data, len(data), replace=True))
        thetaTilde = float(np.median(resamp_data))
        thetaTildeList.append(thetaTilde)
    
    thetaTildeList.sort()
    
    indLower = int(np.round(nresamp * ((1.0 - alpha) / 2))) - 1 # 1 - alpha is, e.g., 95%
    indUpper = nresamp - 1 - indLower - 1
        
    lowerCI = 2 * thetaHat - thetaTildeList[indUpper]
    upperCI = 2 * thetaHat - thetaTildeList[indLower]
        
    return lowerCI, upperCI
