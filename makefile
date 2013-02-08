CC=gcc
COPTS=-O2 -Wall
LINKOPTS=-lm

# Sources
S_HA2=zgbn-ha2.c
S_HA1=zgbn-ha1.c
S_NN=zgb-nn.c

# Ejecutables
E_NN=runNn
E_HA2=runHa2
E_HA1=runHa1

# Archivos miscelaneos para eliminar
MISC_FILES=*.dat

all: nn ha2 ha1

ha1: ${S_HA1}
	${CC} ${COPTS} -o ${E_HA1} ${S_HA1} ${LINKOPTS}

ha2: ${S_HA2}
	${CC} ${COPTS} -o ${E_HA2} ${S_HA2} ${LINKOPTS}

nn: ${S_NN}
	${CC} ${COPTS} -o ${E_NN} ${S_NN} ${LINKOPTS}

clean: clean-dat
	rm -rf ${E_NN} ${E_HA1} ${E_HA2}

clean-dat:
	rm -rf ${MISC_FILES}

cleanbu:
	rm -rf *~
