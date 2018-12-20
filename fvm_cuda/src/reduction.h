
__global__ void finalReduceMin(real *in, real *out ,int n) {
    real val = 1e99;
    for(int i = blockIdx.x*blockDim.x + threadIdx.x; i<n;i+=blockDim.x*gridDim.x) {
        if (in[i] < val)
            val = in[i];
    }
    val = blockReduceMin(val);
    if (threadIdx.x ==0) out[blockIdx.x]=val;
    return;
}
 return val;
}
__global__ void finalReduceSum(real *in, real *out ,int n) {
    real val = 0.;
    for(int i = blockIdx.x*blockDim.x + threadIdx.x; i<n;i+=blockDim.x*gridDim.x) {
    	val += in[i];
    }
    val = blockReduceSum(val);
    if (threadIdx.x ==0) out[blockIdx.x]=val;
    return;
}

__global__ void finalReduceBoolOR(int *in, int *out ,int n) {
    int val = 0;
    for(int i = blockIdx.x*blockDim.x + threadIdx.x; i<n;i+=blockDim.x*gridDim.x) {
    	val |= in[i];
    }
    val = blockReduceBoolOR(val);
    if (threadIdx.x ==0) out[blockIdx.x]=val;
    return;
}





