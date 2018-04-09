CUDA_PATH = /usr/local/cuda-6.5
MYARCH = 20
CC = c++
NVCC = $(CUDA_PATH)/bin/nvcc
CFLAGS = -L$(CUDA_PATH)/lib64 -Wl,-rpath -Wl,$(CUDA_PATH)/lib64 -lcudart -lcublas -O3 -std=c++0x # -fopenmp
NVCCFLAGS= -arch=compute_$(MYARCH) -code=sm_$(MYARCH) -I$(CUDA_SDK_PATH)/C/common/inc -lcublas --ptxas-options=-v -O3

CUSRCS = level3.cu
CSRCS = main_parallel.cpp
EXECNAME = nmf_parallel

CUOBJS = $(CUSRCS:.cu=.o__cu)
COBJS = $(CSRCS:.c=.o__c)

$(EXECNAME): $(CUOBJS) $(COBJS)
	$(CC) $(CFLAGS) $^ -o $@ 

%.o__c: %.c
	$(CC) -o $@ -c $<

%.o__cu: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

clean:
	rm -f core *.o__cu *.o__c *~ $(EXECNAME)
