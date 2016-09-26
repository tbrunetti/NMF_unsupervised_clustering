# NMF_unsupervised_clustering

All the files and scripts in this directory are made to cluster data using NMF and unsupervised learning techniques.  Each module/script is fully functional by itself, however, for convenience and work flow, a BASH script has been provided to streamline all the modules together in one call.

#### Dependencies
------------------
In order to run the entire suite of modules the following dependencies are required:
* numpy (http://docs.scipy.org/doc/numpy-1.10.1/user/install.html)
* matplotlib (http://matplotlib.org/users/installing.html)
* pandas (http://pandas.pydata.org/pandas-docs/stable/install.html)
* seaborn (https://stanford.edu/~mwaskom/software/seaborn/installing.html)

#### Installation and Configuration
------------------------------------
Installation is not necessary, just clone the repository into a directory and it is ready for use.  The streamlined run_NMF.sh may need to be turned into an executable.
```
git clone https://github.com/tbrunetti/NMF_unsupervised_clustering.git
chmod a+x run_NMF.sh
```
If the user chooses to run the streamlined suite of modules, the run_NMF.sh should be configured with the location and specification of the required variables located at the head of the file.  Otherwise, each individual module is equipped with argparse and the -h flag should denote to the user which arguments are required. NOTE: When running the script, you must have access to user X session.


#### Running run_NMF.sh
------------------------

#### Output Files
------------------
