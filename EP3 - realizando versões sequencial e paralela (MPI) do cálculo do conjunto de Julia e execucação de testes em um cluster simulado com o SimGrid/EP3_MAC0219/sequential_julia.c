#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>




/*
 * compute_julia_pixel(): calcula os valores RGB de um pixel em
 *                        uma imagem específica de um conjunto de Julia.
 *
 *  Entrada:
 *      (x, y):           coordenadas do pixel
 *      (largura, altura):  dimensões da imagem
 *      tint_bias:        um valor float para ajustar a tonalidade (1.0 é "sem ajuste")
 *  Saída:
 *      rgb: um array já alocado de 3 bytes onde serão escritos os valores
 *           R, G e B.
 *
 *  Retorno:
 *      0 em caso de sucesso, -1 em caso de falha
 *
 */

int compute_julia_pixel(int x, int y, int largura, int altura, float tint_bias, unsigned char *rgb) {

  // Verifica se as coordenadas são válidas
  if ((x < 0) || (x >= largura) || (y < 0) || (y >= altura)) {
    fprintf(stderr, "Coordenadas inválidas (%d,%d) para um pixel em uma imagem de %d x %d\n", x, y, largura, altura);
    return -1;
  }

  // "Amplia" a visualização para mostrar uma área agradável do conjunto de Julia
  float X_MIN = -1.6, X_MAX = 1.6, Y_MIN = -0.9, Y_MAX = +0.9;
  float float_y = (Y_MAX - Y_MIN) * (float)y / altura + Y_MIN;
  float float_x = (X_MAX - X_MIN) * (float)x / largura + X_MIN;

  // Ponto que define o conjunto de Julia
  float julia_real = -.79;
  float julia_img = .15;

  // Número máximo de iterações
  int max_iter = 300;

  // Calcula a convergência da série complexa
  float real = float_y, img = float_x;
  int num_iter = max_iter;
  while ((img * img + real * real < 2 * 2) && (num_iter > 0)) {
    float xtemp = img * img - real * real + julia_real;
    real = 2 * img * real + julia_img;
    img = xtemp;
    num_iter--;
  }

  // Pinta o pixel com base no número de iterações usando uma coloração estilizada
  float color_bias = (float) num_iter / max_iter;
  rgb[0] = (num_iter == 0 ? 200 : -500.0 * pow(tint_bias, 1.2) * pow(color_bias, 1.6));
  rgb[1] = (num_iter == 0 ? 100 : -255.0 * pow(color_bias, 0.3));
  rgb[2] = (num_iter == 0 ? 100 : 255 - 255.0 * pow(tint_bias, 1.2) * pow(color_bias, 3.0));

  return 0;
}

/* write_bmp_header():
 *
 *   Entrada:
 *      f: Um arquivo aberto para escrita ('w') 
 *      (largura, altura): dimensões da imagem
 *   
 *   Retorno:
 *      0 em caso de sucesso, -1 em caso de falha
 *
 */

int write_bmp_header(FILE *f, int largura, int altura) {

  unsigned int row_size_in_bytes = largura * 3 + 
	  ((largura * 3) % 4 == 0 ? 0 : (4 - (largura * 3) % 4));

  // Define todos os campos no cabeçalho do BMP
  char id[3] = "BM";
  unsigned int filesize = 54 + (int)(row_size_in_bytes * altura * sizeof(char));
  short reserved[2] = {0,0};
  unsigned int offset = 54;

  unsigned int size = 40;
  unsigned short planes = 1;
  unsigned short bits = 24;
  unsigned int compression = 0;
  unsigned int image_size = largura * altura * 3 * sizeof(char);
  int x_res = 0;
  int y_res = 0;
  unsigned int ncolors = 0;
  unsigned int importantcolors = 0;

  // Escreve os bytes no arquivo, mantendo o controle do
  // número de "objetos" escritos
  size_t ret = 0;
  ret += fwrite(id, sizeof(char), 2, f);
  ret += fwrite(&filesize, sizeof(int), 1, f);
  ret += fwrite(reserved, sizeof(short), 2, f);
  ret += fwrite(&offset, sizeof(int), 1, f);
  ret += fwrite(&size, sizeof(int), 1, f);
  ret += fwrite(&largura, sizeof(int), 1, f);
  ret += fwrite(&altura, sizeof(int), 1, f);
  ret += fwrite(&planes, sizeof(short), 1, f);
  ret += fwrite(&bits, sizeof(short), 1, f);
  ret += fwrite(&compression, sizeof(int), 1, f);
  ret += fwrite(&image_size, sizeof(int), 1, f);
  ret += fwrite(&x_res, sizeof(int), 1, f);
  ret += fwrite(&y_res, sizeof(int), 1, f);
  ret += fwrite(&ncolors, sizeof(int), 1, f);
  ret += fwrite(&importantcolors, sizeof(int), 1, f);

  // Sucesso significa que escrevemos 17 "objetos" com êxito
  return (ret != 17);
}
void save_to_file(unsigned char *pixels, int largura, int altura)
{
    FILE *file = fopen("julia.bpm", "w");
    write_bmp_header(file,largura,altura);
    // Escrevendo os pixels após o cabeçalho
    int y;
    int x;
    for (y=0; y < altura; y++) {
        for (x=0; x < largura; x++) {
            fwrite(&(pixels[y * 3 * largura + x * 3]), sizeof(char), 3,file);
        }
        // padding no caso de um número par de pixels por linha
    unsigned char padding[3] = {0,0,0};
    fwrite(padding, sizeof(char), ((largura * 3) % 4),file);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Uso: %s <n> \n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int largura=2*n;
    int altura=n;
    int len = 3*largura*altura;
    unsigned char *vetor = (unsigned char *)malloc(len * sizeof(unsigned char));
    
    int indice;
    clock_t inicio, fim; 
    double tempo_decorrido;

    inicio=clock();
    for (int i = 0; i < altura; i++) {
        for (int j = 0; j < largura; j++) {
            indice = 3 * (i * largura + j);
            compute_julia_pixel(j, i, largura, altura, 1.0, &vetor[indice]);
        }
    }
    fim=clock();

    tempo_decorrido = (double)(fim - inicio) / CLOCKS_PER_SEC;

    printf("Tempo de execução: %f segundos\n", tempo_decorrido);
    save_to_file(vetor, largura, altura);
   

    free(vetor);
   
    return 0;
}