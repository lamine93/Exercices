CFLAGS = -O3 -arch=sm_35  -Wno-deprecated-gpu-targets
CC = nvcc

#insert the ignore for deprecated gpu targets if using CUDA 8.0
CFLAGS += $(shell if [ `nvcc --version | grep 8.0 | wc -l` != 0 ]; then echo -Wno-deprecated-gpu-targets ; fi)

all: laplace2d

laplace: laplace2d.o
	$(CC) $(CFLAGS) -o laplace2d  

clean:
	rm -f *.o laplace2d

%.o: %.cc
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.cu
	$(CC) $(CFLAGS) -o $@ -c $<
