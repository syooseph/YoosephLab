#!/bin/bash

m=`perl -MSet::Scalar -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing Set::Scalar"
	perl -MCPAN -e 'install Set::Scalar'
fi

m=`perl -MGetopt::Long -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing Getopt::Long"
	perl -MCPAN -e 'install Getopt::Long'
fi

m=`perl -MTime::HiRes -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing Time::HiRes"
	perl -MCPAN -e 'install Time::HiRes'
fi

m=`perl -Mconstant -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing constant"
	perl -MCPAN -e 'install constant'
fi

m=`perl -MStatistics::Lite -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing Statistics::Lite"
	perl -MCPAN -e 'install Statistics::Lite'
fi

m=`perl -MPOSIX -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing POSIX"
	perl -MCPAN -e 'install POSIX'
fi

m=`perl -MFile::Basename -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing File::Basename"
	perl -MCPAN -e 'install File::Basename'
fi

m=`perl -MFile::Spec -e 1 2>&1`
if [[ $m != "" ]]; then 
	echo "Installing File::Spec"
	perl -MCPAN -e 'install File::Spec'
fi

# m=`perl -MBio::Lite -e 1 2>&1`
# if [[ $m != "" ]]; then 
# 	echo "Installing Bio::Lite"
# 	perl -MCPAN -e 'install Bio::Lite'
# fi
