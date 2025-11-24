# Test MKL performance
cat('\n=== MKL Performance Test ===\n')
cat('R Version:', R.version.string, '\n')
cat('BLAS:', extSoftVersion()['BLAS'], '\n')
cat('MKL_NUM_THREADS:', Sys.getenv('MKL_NUM_THREADS'), '\n\n')

# Quick benchmark
n <- 3000
A <- matrix(rnorm(n*n), n, n)
cat('Testing 3000x3000 matrix multiplication...\n')
time <- system.time(B <- crossprod(A))['elapsed']
gflops <- (2*n^3/1e9)/time
cat('Time:', time, 'seconds\n')
cat('Performance:', round(gflops, 2), 'GFLOPS\n')

if(time < 2) {
    cat('\n??? MKL appears to be working!\n')
} else {
    cat('\n??? Performance suggests MKL may not be active\n')
}
