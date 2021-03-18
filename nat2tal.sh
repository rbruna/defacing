#!/bin/bash

# Based on FreeSurfer scripts:
# * talairach_avi 1.9 by Nick Schmansky

# Lists the MRI images in the 'input' folder.
IFS=$' \n'
pushd input &> /dev/null
files=( $( ls -1 *.nii.gz 2> /dev/null ) )
stat=$?
popd &> /dev/null

if [ $stat -ne 0 ]
then
  echo "No nii.gz files in the 'input' folder."
  #exit -1
fi

# Creates the working and output directories.
mkdir working &> /dev/null
mkdir output  &> /dev/null

# Defines the environmental variables.
export REFDIR=${FREESURFER_HOME}/average
export MPR2MNI305_TARGET=711-2C_as_mni_average_305

# Operates for each MRI.
for file in ${files[@]}
do
  
  # If no MRI file goes silently to the next file.
  if [ ! -f input/${file} ]
  then
    continue
  fi
  
  # Moves the original MRI file to the 'working' folder.
  mv input/${file} working/
  
  # Enters the working directoy.
  pushd working &> /dev/null
  
  # Extracts the subject name.
  subject=${file/.nii.gz}
  
  
  # Converts the input image to Analyze format.
  echo "Working with subject ${subject}."
  echo "Working with subject ${subject}." 1>> ${subject}.log
  mri_convert ${file} ${subject}.img --conform 1>> ${subject}.log
  
  # Finds the transformation to Talairach space.
  echo "  Transforming the MRI to Talairach space."
  echo "  Transforming the MRI to Talairach space." 1>> ${subject}.log
  mpr2mni305 ${subject} 1>> ${subject}.log
  
  # Saves the transformation matrix.
  echo "  Saving the transformation matrix."
  echo "  Saving the transformation matrix." 1>> ${subject}.log
  avi2talxfm ${subject}.img ${FREESURFER_HOME}/average/mni305.cor.mgz ${subject}_to_711-2C_as_mni_average_305_t4_vox2vox.txt ${subject}_nat2tal.xfm 1>> ${subject}.log
  
  # Moves the output to the output directory.
  mv ${subject}_nat2tal.xfm ../output/
  mv ${subject}.nii.gz ../output/
  mv ${subject}.log ../output/
  
  # Deletes all the working data.
  echo "  Deleting working files."
  rm ${subject}*
  
  # Goes back to the main directory.
  popd &> /dev/null
  
done

# If no working session (i.e. empty working folder) deletes the working folder.
content=$( find working/* ! -type l )
if [[ ! -n ${content[@]} ]]
then
  rm -r working
fi
