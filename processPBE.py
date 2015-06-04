import numpy as np

def moment(C_list_ts, v_list, order):
	out = np.sum(C_list_ts * v_list**order, axis=1)
	return out


def monomer_conservation(M1_ts):
	ini_monomerC = M1_ts[0]
	fin_monomerC = M1_ts[-1]
	conservation = fin_monomerC / ini_monomerC

	results = np.array([ini_monomerC, fin_monomerC, conservation])
	return results

def free_monomers(C_list_ts):
           return C_list_ts[:,0]

def meanV(C_list_ts, v_list):
           M0_ts = moment(C_list_ts, v_list, 0)
           M1_ts = moment(C_list_ts, v_list, 1)
           return M1_ts / M0_ts

def geoMeanV(C_list_ts, v_list):
           numerator = np.sum( C_list_ts * np.log(v_list), axis=1 )
           denumerator = np.sum(C_list_ts, axis=1)
           return np.exp(numerator / denumerator)
           
def geoStdV(C_list_ts, v_list):
  geoMean = geoMeanV(C_list_ts, v_list)
  out = np.empty_like(geoMean)
  for i in np.arange(C_list_ts.shape[0]):
    out[i] = np.sum(C_list_ts[i] * (np.log(v_list) - np.log(geoMean[i]))**2) / (np.sum(C_list_ts[i]) - 1)
  out = np.exp(np.sqrt(out))
  return out
  
           
def meanR(C_list_ts, v_list, R_list):
    numerator   = np.sum(C_list_ts * R_list, axis=1) 
    denominator = moment(C_list_ts, v_list, 0)
    return numerator / denominator

def rmsR(C_list_ts, v_list, R_list):
    numerator = np.sum(C_list_ts * (R_list * v_list)**2, axis=1)
    denominator = moment(C_list_ts, v_list, 2)
    return np.sqrt(numerator / denominator)


def monomer_distribution(C_list_ts, v_list):
    """Calculates the distribution of monomers in cluster."""
    numerator = C_list_ts * v_list
    denominator = moment(C_list_ts, v_list, 1)
    return numerator / denominator

def fitting_stat(x, y):
    slope, intercept, r, prob2, see = stats.linregress(x, y)

    if len(x) > 2:
        see=see*np.sqrt(len(x)/(len(x)-2.))

        mx = x.mean()
        sx2 = ((x-mx)**2).sum()
        sd_intercept = see * np.sqrt(1./len(x) + mx*mx/sx2)
        sd_slope = see * np.sqrt(1./sx2)

    results=np.zeros(5)

    results[0]=slope
    results[1]=intercept
    results[2]=r
    if len(x) > 2:
        results[3]=sd_slope
        results[4]=sd_intercept

    return results

