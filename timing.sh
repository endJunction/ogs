parallel -j1 -u 'git checkout 2dGroundWaterFlow~{} && make -j2 SimpleGroundWaterFlow && ../2dGW_test.sh /data/OGS6_test_data/GroundwaterFlow_Quad4/1e6Cells/case -t > ../test_revs/`git rev-parse HEAD` 2>&1' ::: `seq 32 -1 0`

parallel -j 1 -u ./process_perf_output.sh ::: `git rev-list --reverse --max-count=34 2dGroundWaterFlow` > 1e6Cells.times

plot "1e6Cells.times" u 3 w lp, "1e6Cells.times" u 4 w lp, "1e6Cells.times" u 5 w lp axes x1y2

for i in `seq 1 10`; do ./2dGW_test.sh /data/OGS6_test_data/GroundwaterFlow_Quad4/Quads-x200-y100/case -t 2>&1; done | sed 's/_/ /g' | awk 'BEGIN{N=10}{S=S+$3; U=U+$4; M=M+$6};END{print S/N" " U/N " " M/1024/N}'
#!/bin/bash

for N_MAT_LEN in `seq 0 1000`; do
    echo -n "${N_MAT_LEN} "
    sed -i "s/N_MATERIAL_LENGTH=[0-9]*/N_MATERIAL_LENGTH=$N_MAT_LEN/" SimpleTests/FemTests/CMakeFiles/SimpleGroundWaterFlow.dir/flags.make
    make -j2 SimpleGroundWaterFlow >/dev/null 2>&1
    for i in `seq 1 10`; do
        ./2dGW_test.sh /data/OGS6_test_data/Quads-x200-y100/case -t 2>&1
    done \
        | grep ^sys_user \
        | sed 's/_/ /g' \
        | awk 'BEGIN{N=10}{S=S+$3; U=U+$4; M=M+$6};END{print S/N" " U/N " " M/1024/N}'
done
