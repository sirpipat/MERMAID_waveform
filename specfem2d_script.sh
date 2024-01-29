#!/bin/sh

# Set up and run SPECFEM2D simulation remotely
#
# Usage: specfem2d_script <remotedir> <specfembin>
#
# Make sure that all directories do NOT end with a trailing slash '/'
#
# Last modified by sirawich-at-princeton.edu, 02/10/2023

#####################################################################
# Help                                                              #
#####################################################################
Help()
{
    # Display Help
    echo "Set up and run SPECFEM2D simulation remotely"
    echo
    echo "Usage: specfem2d_script [-d|h] <remotedir> <specfembin>"
    echo
    echo "required arguments":
    echo "remotedir        directory containing DATA/datafiles"
    echo "specfembin       directory to SPECFEM2D binaries"
    echo "Make sure that all directories do NOT end with a trailing slash '/'"
    echo 
    echo "options:"
    echo "d     Debug level"
    echo "      0 - keep only seismograms and source-time function"
    echo "      1 - also keep small information files about the run"
    echo "      2 - also keep JPEG snapshots"
    echo "      3 - keep everything"
    echo "h     Print this help"
}

#####################################################################
# Main Program                                                      #
#####################################################################

# Reading arguments with getopts options
while getopts ":d:h" option; do
    case $option in
	d) # debug level
	    echo Debug level = $OPTARG
	    debug=$OPTARG;;
	h) # display Help
	    Help
	    exit;;
       \?) # invalid option
	    echo "Error: Invalid option"
	    exit;;
    esac
done

# Remove all options passed by getopts options
shift "$(($OPTIND -1))"

# where this function is called
calldir=$(pwd)
echo calling specfem2d_script from $calldir

# move to that remote directory
cd $1

# create a folder for the outputs
echo making OUTPUT_FILES/ directory in $1
mkdir -p OUTPUT_FILES/

# remove any files in OUTPUT_FILES/ if any exist
echo removing files inside ${1}/OUTPUT_FILES/
rm -f OUTPUT_FILES/*

# remove any existing executables
echo removing any existing executables
rm -f xmeshfem2D
rm -f xspecfem2D

# link the executables from the specfembin
echo linking the executables to the original file
if [ ${2:0:1} = "/" ]
then
# Full path starting from root
ln -s ${2}/xmeshfem2D .
ln -s ${2}/xspecfem2D .
else
# Assuming relative path (relative to calling directory)
ln -s ${calldir}/${2}/xmeshfem2D .
ln -s ${calldir}/${2}/xspecfem2D .
fi

# Run xmeshfem2D
echo running xmeshfem2D...
./xmeshfem2D

# Run xmeshfem2D
echo running xspecfem2D...
./xspecfem2D

# Clean up the executable
echo clean up the executables
rm -f xmeshfem2D
rm -f xspecfem2D

# Remove the output files
echo Debug level = $debug
if [ $debug -lt 1 ] ; then
    # keep only seismogram and source-time function
    echo removing the information files 
    rm ${1}/OUTPUT_FILES/*.gnu
    mv ${1}/OUTPUT_FILES/plot_source_time_function.txt ${1}/OUTPUT_FILES/plot_source_time_function.temp
    rm ${1}/OUTPUT_FILES/*.txt
    mv ${1}/OUTPUT_FILES/plot_source_time_function.temp ${1}/OUTPUT_FILES/plot_source_time_function.txt
fi
if [ $debug -lt 2 ] ; then
    # remove JPEG snapshots, keep only seismograms
    echo removing the snapshots
    rm ${1}/OUTPUT_FILES/*.jpg
fi
if [ $debug -lt 3 ] ; then
    # remove the model files
    echo removing the model files
    rm ${1}/DATA/*.bin
    rm ${1}/OUTPUT_FILES/*.bin
    rm ${1}/OUTPUT_FILES/*.vtk
fi

echo The run is complete.
