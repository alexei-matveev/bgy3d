WORK_DIR = 
CMDFLAG = -N 32 -rho 0.018 -beta 1.1989  -pair -norm_tol 1.0e-2 \
	  -max_iter 10000 -L 10.0 -zpad 10.0 -damp_start 1.0 -lambda 0.02
EXE =


BGY2site:
	-@cd ${WORK_DIR}; \
	    ${EXE} -BGY2site ${CMDFLAG}
	-@cd $(CURDIR)

BGYM2site:
	-@cd ${WORK_DIR}; \
	    ${EXE} -BGYM2site ${CMDFLAG}

clean:
	rm -rf ${WORK_DIR}
