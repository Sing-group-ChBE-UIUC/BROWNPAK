#!/usr/bin/env python

#BE CAREFUL WITH THIS ROUTINE.
#
#Syncs fn_traj and fn_stats to make sure that the last record corresponds to
#nts <= nts (in fn_revive).

import argparse
import os
import struct

parser = argparse.ArgumentParser()
parser.add_argument("fr", help="Name of the revive file")
parser.add_argument("--ft", help="Name of the trajectory file")
parser.add_argument("--fs", help="Name of the statistics file")

args = parser.parse_args()

#Name of revive file
assert os.path.exists(args.fr)
fn_revive = args.fr
print('fn_revive = %s'%fn_revive)

#Name of trajectory file.
if args.ft:
    assert os.path.exists(args.fr)
    fn_traj = args.ft
    print('fn_traj = %s'%fn_traj)
else:
    fn_traj = None

#Name of stats file.
if args.fs:
    assert os.path.exists(args.fs)
    fn_stats = args.fs
    print('fn_stats = %s'%fn_stats)
else:
    fn_stats = None

with open(fn_revive, 'rb') as fh_revive:
    fh_revive.read(4) #This is leql (4 bytes), throw away
    buf = fh_revive.read(8) #This is nts, we need to decode this as int64 (long long).
    ubuf = struct.unpack('<q', buf)
    nts_revive = ubuf[0]
    print('nts_revive = %d'%nts_revive)

if fn_traj:
    print('Checking if fn_traj is good.')
    with open(fn_traj, 'rb') as fh_traj:
        #Get header size
        buf = fh_traj.read(4)
        ubuf = struct.unpack('<i', buf)
        hs = ubuf[0] #Header size
        print('fn_traj: Header size (bytes) = %d'%hs)
    
        #Get frame size
        buf = fh_traj.read(4)
        ubuf = struct.unpack('<i', buf)
        fs = ubuf[0] #Frame size
        print('fn_traj: Frame size (bytes) = %d'%fs)
    
        #Number of frames
        num_frames = (os.path.getsize(fn_traj)-hs)//fs
        print('fn_traj: Number of frames = %d'%num_frames)
        if num_frames == 0:
            #No frames, hence no need to truncate
            trunc_traj = False
            print('fn_traj is OK.')
        else:
            #Check if the last frame has nts <= nts_revive
            trunc_traj = True
            offset = hs + (num_frames-1)*fs
            fh_traj.seek(offset, 0)
            buf = fh_traj.read(8)
            ubuf = struct.unpack('<q', buf)
            nts = ubuf[0]
            print('fn_traj: Last frame, nts = %d'%nts)
            if nts <= nts_revive:
                #Last frame has nts <= nts_revive, hence no need to truncate
                trunc_traj = False
                print('fn_traj is OK.')
            else:
                print('fn_traj is bad, will truncate.')
    
    #Truncate fn_traj
    if trunc_traj:
        #Find out offset for truncation. num_frames is > 0.
        fh_traj = open(fn_traj, 'r+b')
        #Loop backwords over the number of frames. The lower bound is 1 in case
        #the file needs to be truncated just after the header.
        iframe = num_frames
        while iframe > 0:
            if (iframe == 1) or (nts <= nts_revive):
                #Truncate the file after current position.
                print('Truncating fn_traj')
                fh_traj.truncate()
                #Close the file. 
                fh_traj.close()
                break
            offset = hs + (iframe-1)*fs
            fh_traj.seek(offset, 0)
            buf = fh_traj.read(8)
            ubuf = struct.unpack('<q', buf)
            nts = ubuf[0]
            iframe -= 1


if fn_stats:
    #Check fn_stats
    print('Checking if fn_stats is good.')
    with open(fn_stats,'r') as fh_stats:
        all_lines = fh_stats.readlines()
    hdr = all_lines[0].rstrip('\n')
    num_hdr_cols = len(hdr.split())
    #Check last line
    records = all_lines[-1].rstrip('\n').split()
    num_cols = len(records)
    if num_cols == num_hdr_cols:
        nts = int(records[0])
        if nts <= nts_revive:
            print('fn_stats OK.')
            create_fs = False
        else:
            print('fn_stats is bad, will chop off.')
            create_fs = True
    else:
        print('fn_stats is bad, will chop off.')
        create_fs = True

    if create_fs: 
        with open(fn_stats,'w') as fh_stats:
            fh_stats.write(hdr+'\n')
        
            for line in all_lines[1:]:
                records = line.rstrip('\n').split()
                num_cols = len(records)
                if num_cols == num_hdr_cols:
                    nts = int(records[0])
                    if nts <= nts_revive:
                        fh_stats.write(line)
