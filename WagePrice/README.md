# Wage Price Spiral, Lorenzoni and Werning (2023)
* wpspiral.nb replicates Figure 2 from the paper
    * Note: y-axis in current version (Brookings, 2023) is mislabeled
* wpmulti.nb gives extension where we consider an isolated supply shock in a multigood economy
    * Assumes expenditure share remains constant and GE effect of shock to mpl doesn't create further price distortions
    * Second assumption critical, still thinking through how to relax this
    * Turns out under these assumptions, fundamentally the same as single good case. Plan to upload derivations at some point
* Have also developed an open economy, multisector extension
    * Extremely complicated even in this vanilla setting
    * Have to [solve](https://springerplus.springeropen.com/articles/10.1186/s40064-016-1810-8#Tab1) a non-linear Fredolm integral of the second kind (if you want to stick with continuous time)
    * Code (python) available upon request          
