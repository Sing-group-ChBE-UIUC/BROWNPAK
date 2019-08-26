#!/usr/bin/env bash

#Driver for program extract.f90.

traj_home="$HOME/workspace/dev/brownpak-0.1/bin"

traj_dir=$traj_home

itraj=""
#if itraj is set to zero-length string
if [ -z $itraj ]; then
    fn_cfg="$traj_dir/chns-10.cfg"
    fn_traj="$traj_dir/traj.bin"
else
    fn_cfg="$traj_dir/lin-30.cfg.$itraj"
    fn_traj="$traj_dir/traj.bin.$itraj"
fi

frame_dir="frmdir"

ifrm_beg=1
ifrm_end=-1
ifrm_stp=1

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
./post "$fn_cfg" "$fn_traj" "$frame_dir" "$ifrm_beg" "$ifrm_end" "$ifrm_stp"
