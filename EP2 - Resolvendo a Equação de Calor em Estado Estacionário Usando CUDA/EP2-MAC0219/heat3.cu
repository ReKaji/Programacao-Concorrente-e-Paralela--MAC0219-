/* 
Este arquivo implementa a Tarefa 3 do EP2 da disciplina MAC0219
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define WALL_TEMP 20.0
#define FIREPLACE_TEMP 100.0
#define BODY_TEMPERATURE 37.0
#define FIREPLACE_START 3
#define FIREPLACE_END 7
#define ROOM_SIZE 10

void initialize(double *h,double *seqh,  int n)
{
    int fireplace_start = (FIREPLACE_START * n) / ROOM_SIZE;
    int fireplace_end = (FIREPLACE_END * n) / ROOM_SIZE;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == 0 || i == n - 1 || j == 0 || j == n - 1)
            {
                h[i * n + j] = (i == n - 1 && j >= fireplace_start && j <= fireplace_end) ? FIREPLACE_TEMP : WALL_TEMP;
                seqh[i * n + j] = (i == n - 1 && j >= fireplace_start && j <= fireplace_end) ? FIREPLACE_TEMP : WALL_TEMP;
            }
            else if((i>=n/2-n/10 && i<=n/2+n/10) && 
                    (j>=n/2-n/10 && j<=n/2+n/10)){
                h[i*n+j]=BODY_TEMPERATURE;
                seqh[i*n+j]=BODY_TEMPERATURE;

            }
            else
            {
                h[i * n +j] = 0.0;
                seqh[i * n +j] = 0.0;

            }
        }
    }
}


bool compara_cpu_gpu(double *h, double *seqh, int n){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (fabs(h[i*n+j]-seqh[i*n+j])>0.5) 
                return false;
              
        }
        
    }
    return true;
}

__global__ void jacobi_iteration(double *h, double *g, int n) {
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    int stride_x = gridDim.x * blockDim.x;
    int stride_y = gridDim.y * blockDim.y; 

    for (int row = i; row < n - 1; row += stride_x) {
        for (int col = j; col < n - 1; col += stride_y) {
 
            if (row > 0 && row < n -1 && col > 0 && col < n-1 && (!(row >= n / 2 - n / 10 && row <= n / 2 + n / 10 && 
                      col >= n / 2 - n / 10 && col <= n / 2 + n / 10))) {
                g[row * n + col] = 0.25 * (h[(row - 1) * n + col] +
                                           h[(row + 1) * n + col] +
                                           h[row * n + col - 1] +
                                           h[row * n + col + 1]);
            }
        }
    }
    __syncthreads();
    for (int row = i; row < n -1; row += stride_x) {
        for (int col = j; col < n -1; col += stride_y) {
            if (row > 0 && row < n -1 && col > 0 && col < n-1 && (!(row >= n / 2 - n / 10 && row <= n / 2 + n / 10 && 
                      col >= n / 2 - n / 10 && col <= n / 2 + n / 10)))
                h[row*n+col]=g[row*n+col];
       }}
    __syncthreads();
}

void seq_jacobi_iteration(double *h, double *g, int n, int iter_limit)
{
    for (int iter = 0; iter < iter_limit; iter++)
    {
        for (int i = 1; i < n - 1; i++)
        {
            for (int j = 1; j < n - 1; j++)
            {
                if (i > 0 && i < n -1 && j > 0 && j < n-1 && (!(i >= n / 2 - n / 10 && i <= n / 2 + n / 10 && 
                      j >= n / 2 - n / 10 && j <= n / 2 + n / 10)))
                    g[i * n +j] = 0.25 * (h[(i-1)*n+ j] + h[(i + 1)*n+j] + h[i * n + j - 1] + h[i * n + j + 1]);
            }
        }
        for (int i = 1; i < n - 1; i++)
        {
            for (int j = 1; j < n - 1; j++)
            {
                if (i > 0 && i < n -1 && j > 0 && j < n-1 && (!(i >= n / 2 - n / 10 && i <= n / 2 + n / 10 && 
                      j >= n / 2 - n / 10 && j <= n / 2 + n / 10)))
                    h[i * n + j] = g[i*n+j];
            }
        }
    }
}
double calculate_elapsed_time(struct timespec start, struct timespec end)
{
    double start_sec = (double)start.tv_sec * 1e9 + (double)start.tv_nsec;
    double end_sec = (double)end.tv_sec * 1e9 + (double)end.tv_nsec;
    return (end_sec - start_sec) / 1e9;
}

void save_to_file(double *h, int n)
{
    FILE *file = fopen("room.txt", "w");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(file, "%lf ", h[i*n+j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        fprintf(stderr, "Uso: %s <número de pontos> <limite de iterações> <t> <b>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int iter_limit = atoi(argv[2]);
    int t= atoi(argv[3]);
    int b=atoi (argv[4]);

    double *h = (double *)malloc(n * n * sizeof(double));
    double *g = (double *)malloc(n * n * sizeof(double));

    double *seqh = (double *)malloc(n * n * sizeof(double));
    double *seqg = (double *)malloc(n * n * sizeof(double));
    

    double *dh, *dg;

    
    cudaMalloc((void**)&dh, sizeof(double) * n * n);
    cudaMalloc((void**)&dg, sizeof(double) * n * n);
    
    



    if (h == NULL || g == NULL)
    {
        fprintf(stderr, "Erro ao alocar memória para h ou g\n");
        exit(EXIT_FAILURE);
    }


   
    initialize(h, seqh,n);

    float tempo_host_device=0;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    cudaMemcpy (dg, g,n*n*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy (dh, h, n*n*sizeof(double), cudaMemcpyHostToDevice);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&tempo_host_device, start, stop);
    tempo_host_device/=1000;
    printf("Tempo de cudaMemcpyHostToDevice de g e h : %.9f segundos\n", tempo_host_device);

    int raiz_t= (int)sqrt(t);
    int raiz_b=(int)sqrt(b);
    int tam_bloco1;
    int tam_bloco2;
    int tam_thread1;
    int tam_thread2;

    if (raiz_t * raiz_t == t){
        tam_thread1=  raiz_t;
        tam_thread2= raiz_t;
    }
    else
        {
            t=t/2;
            int raiz= (int)sqrt(t);
            tam_thread1= raiz;
            tam_thread2=2* raiz;
        }

    if (raiz_b * raiz_b == b){
        tam_bloco1=  raiz_b;
        tam_bloco2= raiz_b;
    }
    else
        {
            b=b/2;
            int raiz_b= (int)sqrt(b);
            tam_bloco1= raiz_b;
            tam_bloco2=2*raiz_b;
        }




    dim3 threads (tam_thread1,tam_thread2);
    dim3 blocks(tam_bloco1,tam_bloco2);

    float tempo_cuda = 0;
    float tempo_tarefa=0;
    for (int i=0; i< iter_limit;i++){
        
        cudaEventRecord(start);
        jacobi_iteration<<<blocks,threads>>>(dh, dg, n);
        cudaDeviceSynchronize();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&tempo_tarefa, start, stop);
        tempo_cuda+=tempo_tarefa;
    }
    tempo_cuda/=1000;
    printf("Tempo de execução versão CUDA: %.9f segundos\n", tempo_cuda);
    
    float tempo_device_host=0;

    cudaEventRecord(start);
    cudaMemcpy (h, dh, n*n*sizeof(double), cudaMemcpyDeviceToHost);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&tempo_device_host, start, stop);
    tempo_device_host/=1000;
    printf("Tempo de cudaMemcpyDeviceToHost de h: %.9f segundos\n", tempo_device_host);
   
    struct timespec start1, end1;

    clock_gettime(CLOCK_MONOTONIC, &start1);
    seq_jacobi_iteration(seqh, seqg, n, iter_limit);
    clock_gettime(CLOCK_MONOTONIC, &end1);

    save_to_file(h, n);

    double elapsed_time = calculate_elapsed_time(start1, end1);
    printf("Tempo de execução versão sequncial: %.9f segundos\n", elapsed_time);
    if (compara_cpu_gpu(h, seqh, n)) 
        printf("As versões da GPU e CPU produzem o mesmo resultado\n");
    else 
        printf("As versões da GPU e CPU NÃO produzem o mesmo resultado\n");


    free(h);   
    free(g);
    free(seqg);
    free(seqh);
    cudaFree(dg);
    cudaFree(dh);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return 0;
}