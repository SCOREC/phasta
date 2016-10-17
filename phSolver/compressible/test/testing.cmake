set(casename staticBubble)
set(CDIR ${CASES}/${casename}/run)
set(NPROCS 8)
c_serial_test(inpCfg_${casename} ${CDIR} cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
c_serial_test(resCfg_${casename} ${CDIR} cp -r ${CDIR}/../chef/${NPROCS}-procs_case ${CDIR})
c_parallel_test(staticBubble ${NPROCS} ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
