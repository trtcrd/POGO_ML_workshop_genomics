#! /bin/sh

# adding the last sources for R 
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
# For R core
sudo apt install -y r-base
# For dependency of devtools in R
sudo apt install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
# For R Studio Server
sudo apt install -y gdebi-core
wget https://download2.rstudio.org/rstudio-server-1.1.383-amd64.deb
sudo gdebi -n rstudio-server-1.1.383-amd64.deb
sudo rm rstudio-server-1.1.383-amd64.deb
# Make Accessible to /usr/local/lib/R/site-library
sudo chmod 777 -R /usr/local/lib/R/site-library

# Change port from 8787 to 80 for user friendly
sudo echo "www-port=80" >> /etc/rstudio/rserver.conf
sudo rstudio-server restart

# install R libs
R -e 'install.packages("vegan", repos="https://stat.ethz.ch/CRAN/")'
R -e 'install.packages("ranger", repos="https://stat.ethz.ch/CRAN/")'
R -e 'install.packages("doMC", repos="https://stat.ethz.ch/CRAN/")'
R -e 'install.packages("BBI", repos="https://stat.ethz.ch/CRAN/")'
R -e 'install.packages("irr", repos="https://stat.ethz.ch/CRAN/")'
R -e 'install.packages("matrixStats", repos="https://stat.ethz.ch/CRAN/")'
