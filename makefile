CC = nvcc
CC_C = gcc
ARCH = sm_70
#CFLAGS = -O3 
GPUFLAGS = -DUSE_GPU=1 -arch=$(ARCH) -fmad=false
CFLAGS = -O0 -g -G
#CFLAGS = -O3  -pg
#CFLAGS = -O0  -g 
#CFLAGS = -Xcompiler "-fsanitize=address -fno-omit-frame-pointer" -O0  -g  
CONLYFLAGS = -O3  
LIBS = -lm 
#OBJS =  main.o  cpu_dem.o device_dem.o  ParticleSystem.o output.o BoundingBox.o solver_output.o TriangleMesh.o BVH.o cJSON.o settings_loader.o

OBJS =  main.o ParticleSystem.o cJSON.o TriangleMesh.o settings_loader.o solver_output.o  BVH.o BoundingBox.o cpu_dem.o device_dem.o

PROGRAM = myDEM3d

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(GPUFLAGS) $(OBJS) $(LIBS)  -o $(PROGRAM)

%.o : %.cu
	$(CC) $(CFLAGS) $(GPUFLAGS) -c $<

%.o : %.c
	$(CC_C) $(CONLYFLAGS) -c $<

clean:
	rm -f $(PROGRAM)  $(OBJS)

