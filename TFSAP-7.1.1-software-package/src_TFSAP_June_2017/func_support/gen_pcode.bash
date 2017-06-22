#!/bin/bash
#
# run matlab in batch mode (silent) to create .p files
#
# use: gen_pcode.sh MATLAB_BIN <files>
#
MATLAB_BIN=$1

for arg in "$@"
do
    ARGS_IN[index]="$arg"
    let "index+=1"
done 

B_ARRAY="${ARGS_IN[@]:1}"
echo "B_ARRAY = ${B_ARRAY}"



echo 'Starting MATLAB process using pcode() :'

unset DISPLAY
${MATLAB_BIN} > matlab_pcode.out 2>&1 <<EOF
pcode ${B_ARRAY};
exit;
EOF

echo 'Finished MATLAB process'

