CC = nvcc
CFLAGS = -O3  -DUSE_GPU=0 -arch=sm_86 -DUSE_GPU=0
LIBS = -lm 
OBJS =  main.o  cpu_dem.o device_dem.o  ParticleSystem.o output.o BoundingBox.o solver_output.o 

PROGRAM = myDEM3d

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(PROGRAM)

%.o : %.cu
	${CC} ${CFLAGS} -c $<

%.o : %.c
	${CC} ${CFLAGS} -c $<

clean:
	rm -f $(PROGRAM)  $(OBJS)

