CC=gcc
COPTS=-O2 -Wall
LINKOPTS=-lm

# Sources
S_HA=zgbn-ha.c
S_HA1=zgbn-ha1.c
S_NN=zgb-nn.c

# Ejecutables
E_NN=runNn
E_HA=runHa
E_HA1=runHa1

# Archivos miscelaneos para eliminar
MISC_FILES=*.dat

all: nn ha ha1

ha: ${S_HA}
	${CC} ${COPTS} -o ${E_HA} ${S_HA} ${LINKOPTS}

ha1: ${S_HA1}
	${CC} ${COPTS} -o ${E_HA1} ${S_HA1} ${LINKOPTS}

nn: ${S_NN}
	${CC} ${COPTS} -o ${E_NN} ${S_NN} ${LINKOPTS}

clean: clean-dat
	rm -rf ${E_NN} ${E_HA} ${E_HA1}

clean-dat:
	rm -rf ${MISC_FILES}
