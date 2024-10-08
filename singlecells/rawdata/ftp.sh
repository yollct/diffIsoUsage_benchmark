#!/bin/bash
# Run this file in bash with this command:  ./filename
HOST=ftp.ebi.ac.uk
USER=anonymous
ftp -pinv $HOST <<EOF
user $USER
cd biostudies/fire/E-MTAB-/678/E-MTAB-6678/Files
binary
mget "meta_ss2.txt"
mget "raw_data_ss2.txt"
disconnect
bye
EOF