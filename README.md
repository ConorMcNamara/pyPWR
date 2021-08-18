# pyPWR

Despite the proliferation of A/B testing in Data Science, Python lacks a singular package to calculate the power of many different types of statistical tests. Conversely, R has the [pwr](https://github.com/heliosdrm/pwr/tree/master/R) library whose sole purpose is to provide calculations for power, effect size, sample size and significance levels. 

The goal of this package is to reverse engineer the R library to provide a Python equivalent of `pwr`. Included are all equivalent functions in the R library, with some changes to make it more Pythonic (e.g., classes instead of inner functions; underscores for names rather than periods), as well as unit tests to ensure that the results are in-line with what the `pwr` library reports. 


## Quick Example
``` 
from pwr_tests import *
power = pwr_2p_test(h=0.3, n=200, sig_level=0.05, alternative='greater')
power['power']
0.9123145
```

Note that similar to the R library, one of `effect_size`, `sample_size`, `sig_level`, or `power` must be left as None. 

## Notes
Due to the fact that `pwr` uses R's [uniroot](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/uniroot) for root solving whereas I used Scipy's [brentq](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html) or [bisect](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html), there are going to be some minor differences in terms of the values reported for `effect_size` or `sig_level`; however, they should be within error of each other (1e-03)
