#!/usr/bin/env bash

#Driver for program sync_files.py.

traj_home="$HOME/workspace/dev/brownpak-0.1/bin"

traj_dir=$traj_home

itraj=""
#if itraj is set to zero-length string
if [ -z $itraj ]; then
    #Name of trajectory file.
    fn_traj="$traj_dir/traj.bin"
    #Name of stats file.
    fn_stats="$traj_dir/stats.txt"
    #Name of revive file
    fn_revive="$traj_dir/revive.bin"
else
    fn_traj="$traj_dir/traj.bin.$itraj"
    fn_stats="$traj_dir/stats.txt.$itraj"
    fn_revive="$traj_dir/revive.bin.$itraj"
fi

#Run
python sync_files.py "$fn_revive" --ft "$fn_traj" --fs "$fn_stats"
