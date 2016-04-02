# g++4.8.1
if [ "$CXX" == "g++" ]; then sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test; fi

# clang 3.4
if [ "$CXX" == "clang++" ]; then sudo add-apt-repository -y ppa:h-rayflood/llvm; fi

sudo apt-get update -qq

# to prevent IPv6 being used for APT
sudo bash -c "echo 'Acquire::ForceIPv4 \"true\";' > /etc/apt/apt.conf.d/99force-ipv4"
# The ultimate one-liner setup for NeuroDebian repository
bash <(wget -q -O- http://neuro.debian.net/_files/neurodebian-travis.sh)
# But we actually want -devel repository just to get matio backport
sudo sed -ie 's,neuro.debian.net/debian ,neuro.debian.net/debian-devel ,g' /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-get update
# Just to get information about available versions
apt-cache policy libmatio-dev
