#!/bin/bash


usage="$(basename "$0") converts binary surface files of a subject to stl format. The converted files is stored in the folder <mesh> in the subject directory
where:
    -h  show this help text
    -s  subject name"

while getopts 'hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;; 
    s)
      subjid=$OPTARG
      mkdir -p $SUBJECTS_DIR/${subjid}/mesh
      # mris_convert is a Freesurfer command that converts surface files to other formats, such as stl and vtk.
      mris_convert $SUBJECTS_DIR/${subjid}/surf/lh.pial $SUBJECTS_DIR/${subjid}/mesh/lh-pial.stl 
      mris_convert $SUBJECTS_DIR/${subjid}/surf/rh.pial $SUBJECTS_DIR/${subjid}/mesh/rh-pial.stl 
      mris_convert $SUBJECTS_DIR/${subjid}/surf/lh.white $SUBJECTS_DIR/${subjid}/mesh/lh-white.stl 
      mris_convert $SUBJECTS_DIR/${subjid}/surf/rh.white $SUBJECTS_DIR/${subjid}/mesh/rh-white.stl 
     ;;
  esac
done


