
# Programação Concorrente e Paralela (MAC0219) — Projetos e competências técnicas

Este repositório reúne os exercícios práticos da disciplina MAC0219 (Universidade de São Paulo). Cada pasta contém uma implementação experimental que explora modelos diferentes de paralelismo e técnicas de otimização com foco em desempenho, reprodutibilidade e análise experimental.

## Estrutura

- `EP1 - Cálculo do Conjunto de Mandelbrot utilizando Pthreads e OpenMP/EP1-MAC0219/`
	- `mandelbrot_seq.c` — implementação sequencial em C
	- `mandelbrot_pth.c` — paralelização com Pthreads
	- `mandelbrot_omp.c` — paralelização com OpenMP
	- `Planilha_EP1_MAC0219.csv` — medições experimentais e dados de benchmark

- `EP2 - Resolvendo a Equação de Calor em Estado Estacionário Usando CUDA/EP2-MAC0219/`
	- `heat2.cu`, `heat3.cu` — implementações orientadas a GPU em CUDA
	- `Makefile` — regras para compilar e executar as experiências

- `EP3 - realizando versões sequencial e paralela (MPI) do cálculo do conjunto de Julia/EP3_MAC0219/`
	- `sequential_julia.c` — versão sequencial em C
	- `1D_parallel_julia.c` — versão distribuída usando MPI
	- `Makefile`

## O que cada projeto faz (resumo técnico)

EP1 — Mandelbrot (C, Pthreads, OpenMP)
- Implementa a geração do conjunto de Mandelbrot em três variantes (sequencial, multithread com Pthreads e paralelismo de alto nível com OpenMP). O trabalho foca em decomposição de carga, balanceamento e minimização de overhead de sincronização.
- Habilidades práticas: design de partição de trabalho, controle de concorrência, otimizações de loop e coleta de métricas reprodutíveis para comparar estratégias.

EP2 — Equação de calor (CUDA)
- Resolve numericamente a equação de calor em estado estacionário aproveitando a arquitetura GPU. As versões exploram organização de threads, uso de memória compartilhada e trade-offs entre paralelismo e custo de comunicação memória/compute.
- Habilidades práticas: escrita e otimização de kernels CUDA, gerenciamento explícito de memória e perfilamento de desempenho para identificar gargalos.

EP3 — Conjunto de Julia (MPI)
- Paraleliza a geração do conjunto de Julia em um modelo de memória distribuída usando MPI. A versão paralela ilustra particionamento de domínio e comunicação entre processos para consolidar resultados.
- Habilidades práticas: design de comunicação (point-to-point e coletivas), cuidado com granularidade de tarefas e análise do overhead de comunicação.

## Como compilar e executar (exemplos)

Observação: os comandos abaixo assumem um sistema Linux com ferramentas de desenvolvimento instaladas (`gcc`, `nvcc` para CUDA, `mpicc`/`mpirun` para MPI). Ajuste conforme seu ambiente.

EP1 — Mandelbrot:
```bash
cd "EP1 - Cálculo do Conjunto de Mandelbrot utilizando Pthreads e OpenMP/EP1-MAC0219"
gcc -O3 mandelbrot_seq.c -o mandelbrot_seq -lm
gcc -O3 -fopenmp mandelbrot_omp.c -o mandelbrot_omp -lm
gcc -O3 mandelbrot_pth.c -o mandelbrot_pth -pthread -lm
# exemplo de execução (parâmetros dependem da implementação)
./mandelbrot_omp 1024 1024 1000
```

EP2 — Equação de calor (CUDA):
```bash
cd "EP2 - Resolvendo a Equação de Calor em Estado Estacionário Usando CUDA/EP2-MAC0219"
make
./heat2
```

EP3 — Julia (MPI):
```bash
cd "EP3 - realizando versões sequencial e paralela (MPI) do cálculo do conjunto de Julia/EP3_MAC0219"
make
./sequential_julia
mpirun -np 4 ./1D_parallel_julia
```



