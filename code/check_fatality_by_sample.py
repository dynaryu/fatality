# read USGS Expo Cat info
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

def compute_probability(total_fatalities):
    """Pager method compute probaility of fatality in each magnitude bin.

    [0,10), [10,10^2), [10^2,10^3), [10^3, 10^4), [10^4, 10^5), 
    [10^5,)

    :param total_fatalities: List of total fatalities in each MMI class.
    :type total_fatalities: int, float

    :returns: Probability of fatality magnitude bin from
        lognorm.cdf(bin, shape=Zeta, scale=total_fatalities)
    :rtype: list(float) """

    nsamples = float(len(total_fatalities))
    cprob = np.ones(len(magnitude_bin)+1)*100.0
    for j, val in enumerate(magnitude_bin):
        cprob[j] = np.sum(total_fatalities < val)/nsamples*100.0
    return np.hstack((cprob[0], np.diff(cprob)))

filename = '/Users/hyeuk/Project/fatality/data/EXPO_CAT_2007_12.csv'
expo_cat = pd.read_csv(filename)

irow = (expo_cat["ISO_code"] == 'ID') & (~(np.isnan(expo_cat.icol(16)))) & (
    expo_cat["eqID"] != 200503281609) & (expo_cat["eqID"] != 200605262253)
# remove two events to avoid any duplication

sel_dat = expo_cat.ix[irow]

# observed fatality 
observed = sel_dat["PAGER_prefShakingDeaths"].values

mag_observed = np.floor(np.log10(observed)) + 1
mag_observed[np.isinf(mag_observed)] = 1 # 0 range

eqID = sel_dat["eqID"].values
mag = sel_dat["PAGER_prefMag"].values
magType = sel_dat["PAGER_prefMagType"].values

# exposed population including rural and urban
# U090: number of people exposed to MMI 9.0 +/- 0.25, or MMI 8.75 to MMI 9.25
mmi_list = range(4, 11)
pop = np.zeros((len(sel_dat), len(mmi_list)))
for i, val in enumerate(mmi_list):
    str_R = 'R%03d' %(val*10)
    str_U = 'U%03d' %(val*10)
    pop[:,i] = (sel_dat[str_R]+sel_dat[str_U]).values

# check threee sampled fatality rate
data_path = '/Users/hyeuk/Project/inasafe/safe/impact_functions/earthquake/itb_bayesian_earthquake_fatality_model_trial'
rate_org = np.load(os.path.join(data_path, 'worden_berngamma_log_fat_rate_inasafe.npy'))
rate_10 = np.load(os.path.join(data_path, 'worden_berngamma_log_fat_rate_inasafe_10.npy'))
rate_50 = np.load(os.path.join(data_path, 'worden_berngamma_log_fat_rate_inasafe_50.npy'))

res_org = np.dot(pop, rate_org.T)
res_10 = np.dot(pop, rate_10.T)
res_50 = np.dot(pop, rate_50.T)

# compute prob
magnitude_bin = np.power(10, range(1, 6), dtype=float)
prob_org = np.zeros((len(sel_dat), len(magnitude_bin)+1))
prob_10 = np.zeros((len(sel_dat), len(magnitude_bin)+1))
prob_50 = np.zeros((len(sel_dat), len(magnitude_bin)+1))
for i in range(len(sel_dat)):
    prob_org[i, :] = compute_probability(res_org[i])
    prob_10[i, :] = compute_probability(res_10[i])
    prob_50[i, :] = compute_probability(res_50[i])

# magnitude
idx_max = prob_org.argmax(axis=1) +1 # 0-1 goes to 1 instead of zero

# bin 5-5.5
bins = np.arange(4.75, 9.5, 0.5)
digitized = np.digitize(mag, bins)
for i in range(1, len(bins)):
    tf = (digitized == i)
    total_freq = np.sum(tf)
    no_match = np.sum(idx_max[tf]==mag_observed[tf])
    no_non_match = np.sum(~(idx_max[tf]==mag_observed[tf]))
    print "%.2f-%.2f: %d out %d matched, %d missed" %(bins[i-1], bins[i], no_match, total_freq, no_non_match)


for i in range(len(sel_dat)):
    plt.figure()
    temp = np.concatenate((prob_org[i,:], prob_10[i,:], prob_50[i,:])).reshape(3,6).T
    plt.stem(range(6), temp[:,0])
    #plt.stem(np.arange(6)+0.2, temp[:,1], linefmt='r-', markerfmt='ro')
    #plt.stem(np.arange(6)+0.4, temp[:,2], linefmt='g-', markerfmt='go')
    plt.grid(1)
    plt.xticks(range(6),("0-10","10-100","10^2-10^3","10^3-10^4","10^4-10^5","10^5+"))
    plt.title('%s, %s%.2f: %d' %(eqID[i], magType[i], mag[i], observed[i]))
    plt.savefig('../plot/check_mag/check_%s.png' %i)
    plt.close()

# looks like 10 will provide about the same result as the original