#!/usr/bin/env bash

#Driver for program extract.f90.

traj_home="$HOME/workspace/dev/brownpak-0.1/bin"

traj_dir=$traj_home

itraj=""
#if itraj is set to zero-length string
if [ -z $itraj ]; then
    #Name of trajectory file.
    fn_traj="$traj_dir/traj.bin"
    #Name of .cfg file specifying the configuration. The atom positions will
    #be overwritten by the data in fn_traj.
    fn_cfg="$traj_dir/sample-rlx.cfg"
else
    fn_traj="$traj_dir/traj.bin.$itraj"
    fn_cfg="$traj_dir/foo.cfg.$itraj"
fi

#Where to output the frames
frame_dir="frmdir"

#Number of the first frame
ifrm_beg=1
#Number of the last frame. -1 indicates the last available frame.
ifrm_end=-1
#Step over how many frames? 1 indicates consecutive frames.
ifrm_stp=1
#Flag for MPCD data
mflag=3
#Write MPCD data?
with_mpcd_atoms=T

#If frame_dir exists, clear frame_dir of all preexisting frames
#If frame_dir does not exist, create frame_dir
if [[ -d "$frame_dir" ]]; then
    rm -rf "$frame_dir"/frame*
    printf "%s emptied \n" "$frame_dir"
else
    mkdir -p "$frame_dir"
    printf "%s created \n" "$frame_dir"
fi

#Run
./extract "$fn_cfg" "$fn_traj" "$frame_dir" "$ifrm_beg" "$ifrm_end" "$ifrm_stp" \
    "$mflag" "$with_mpcd_atoms"
