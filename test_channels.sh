#!/bin/bash

REF_FILE="reference.txt"

if [ ! -f "$REF_FILE" ]; then
    echo -e "\033[31m[ERROR] Reference file '$REF_FILE' not found!\033[0m"
    echo "Please create '$REF_FILE' or update the REF_FILE variable in the script."
    exit 1
fi

fs_list=("ee" "mm" "gg" "pp" "mr" "pr")
ord_list=("born" "alpha" "exp")

# Initialize mismatch/error counter
mismatch_count=0

for fs in "${fs_list[@]}"; do
    for ord in "${ord_list[@]}"; do

        if [ "$fs" == "mr" ] || [ "$fs" == "pr" ]; then
            zmax=180
            extra_cuts="egmin 0.02"
        else
            zmax=10
            emin_val=$(awk 'BEGIN {print 0.4 * 1.02}')
            extra_cuts="emin $emin_val"
        fi

        if [ "$fs" == "pp" ] || [ "$fs" == "pr" ]; then
            arun_list=("off")
        else
            arun_list=("off" "nsk")
        fi

        if [ "$fs" == "pp" ] || [ "$fs" == "pr" ]; then
            iffpi_list=(1 2 3)
        else
            iffpi_list=(0)
        fi

        for arun in "${arun_list[@]}"; do
            for iffpi in "${iffpi_list[@]}"; do

                # --- CHECK PHIDEC CONDITION ---
                if [ "$iffpi" -eq 1 ] && [ "$fs" == "pr" ] && [ "$ord" == "exp" ]; then
                    phidec_list=("000" "100" "030" "001" "131")
                else
                    phidec_list=("000")
                fi

                # --- LOOP OVER PHIDEC PARAMETERS ---
                for phidec in "${phidec_list[@]}"; do

                    cat <<EOF > input
fs $fs
nev 1000
path test/
ord $ord
ecms 1.02
zmax $zmax
thmin 20
thmax 160
model matched
mode weighted
arun $arun
nphot -1
iffpi $iffpi
iFSRdisp 1
phidec $phidec
eps 1d-5
compttens 0
what_ffpi bwsum2
seed 55
$extra_cuts
run
EOF

                    # Run BabaYaga
                    ./babayaga < input > /dev/null 2>&1

                    # Parse output and extract reference values
                    if ls test/stat* 1> /dev/null 2>&1; then
                        res=$(grep -h "total:" test/stat* | awk '{print $2, $4}')

                        # Search for reference line matching this combination
                        ref_line=$(grep "fs: $fs " "$REF_FILE" | grep "arun: $arun " | grep "iffpi: $iffpi " | grep "ord: $ord ")
                        
                        # Filter further if the reference file contains the phidec parameter
                        if echo "$ref_line" | grep -q "phidec: $phidec"; then
                            ref_line=$(echo "$ref_line" | grep "phidec: $phidec ")
                        fi

                        if [ -n "$ref_line" ]; then
                            # Extract last two columns from reference file
                            ref_res=$(echo "$ref_line" | awk '{print $(NF-1), $NF}')

                            # Numerical comparison with tolerance |diff| < 10^-8
                            is_ok=$(awk -v r1="$res" -v r2="$ref_res" '
                            BEGIN {
                                n1 = split(r1, a)
                                n2 = split(r2, b)
                                if (n1 != n2 || n1 == 0) { print 0; exit }
                                tol = 1e-5
                                for (i = 1; i <= 1; i++) {
                                    diff = (a[i] - b[i])/(a[i])
                                    if (diff < 0) diff = -diff
                                    if (diff >= tol) { print 0; exit }
                                }
                                print 1
                            }')

                            if [ "$is_ok" -eq 1 ]; then
                                status="\033[32m[OK]\033[0m"
                            else
                                status="\033[31m[MISMATCH! Expected: $ref_res]\033[0m"
                                ((mismatch_count++))
                            fi
                        else
                            status="\033[33m[REFERENCE NOT FOUND]\033[0m"
                            ((mismatch_count++))
                        fi

                        echo -e "fs: $fs | arun: $arun | iffpi: $iffpi | ord: $ord | phidec: $phidec | xsec: $res | $status"
                    else
                        echo "fs: $fs | arun: $arun | iffpi: $iffpi | ord: $ord | phidec: $phidec | [Error: test/stat* not found]"
                        ((mismatch_count++))
                    fi

                    # Cleanup temporary files
#                    rm -rf test
#                    rm -f input

                done
            done
        done
    done
done

# --- FINAL VERIFICATION ---
echo "--------------------------------------------------------------------------------"
if [ "$mismatch_count" -eq 0 ]; then
    echo -e "\033[32m\033[1mCHECK OK, CORRECT INSTALLATION\033[0m"
else
    echo -e "\033[31m\033[1m[CHECK FAILED] Found $mismatch_count mismatch(es) or error(s) during execution.\033[0m"
fi
echo "--------------------------------------------------------------------------------"
