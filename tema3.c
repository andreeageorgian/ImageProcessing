#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define lineSize 1024

//Structura in care memoram latimea, inaltimea
//, maximum_gray, precum si matricea imaginii.
struct PGMImg{
	int width;
	int height;
	int maximum_gray;
	int **matrix;
};

//Structura in care memoram matricea 
//filtrului, factorul si deplasamentul.
struct Filter{
	int **matrix;
	int factor;
	int displacement;
};

//Functie de alocare de memorie pentru o matrice.
int **allocate_matrix(int width, int height){
	int i;
	int **mat = (int **)calloc(sizeof(int*), height);
	
	for(i = 0; i < height; i++) {
		mat[i] = (int *)calloc(sizeof(int), width);
	}
	
	return mat;
}

//Functie de realocare de memorie pentru o matrice.
void reallocate_matrix(int **matrix, int height){
	int i;
	
	matrix = (int **)realloc(matrix, height * sizeof(int*));
}

//Functie de eliberare de memorie pentru o matrice.
void deallocate_matrix(int **matrix, int height){
    int i;

    for (i = 0; i < height; ++i) {
        free(matrix[i]);
    }

    free(matrix);
}

//Functie ce citeste matricea, factorul si deplasamentul
//filtrului din fisierul corespunzator.
struct Filter readFilter(const char *file_name){
	struct Filter filter;
	FILE *filter_file = fopen(file_name, "r");
	
	int i = 0 , j = 0, aux;
	char line[20];
	char *token;
	const char s[] = " \n/+";

	if(filter_file == NULL){
		printf("Eroare citire filtru!\n");
	}

	filter.matrix = allocate_matrix(4,4);

	while(fgets(line,sizeof line, filter_file) != NULL){
		token = strtok(line, s);
		filter.matrix[i][j] = atoi(token);
		while( token != NULL ) {
			j++;
      		token = strtok(NULL, s);
      		if(j < 3){
      			filter.matrix[i][j] = atoi(token);
      		} else if(j == 3 && token != NULL){
      			filter.factor = atoi(token);
      		} else if(j == 4 && token != NULL){
      			filter.displacement = atoi(token);
      		}
   		}
   		i++;
   		j = 0;
	}

	fclose(filter_file);
	return filter;
}

//Functie ce completeaza vectorul de vecini al fiecarui nod.
void set_neighbours_array(int* neighbours, char* line, int array_size){
	int i, proc;
	const char s[] = " :\n";
	char *token;

	token = strtok(line, s);
	proc = atoi(token);
	while(token!=NULL){
		if(atoi(token) != proc){
			neighbours[atoi(token)] = 1;
		}
		token = strtok(NULL, s);		
	}
}

//Functie ce extrage linia cu numarul corespunzator
//din fisierul primit ca parametru.
char* get_line(const char* file_in, int line_number) {
	char* line = (char *)calloc(lineSize, sizeof(char));
	FILE *file = fopen(file_in, "r");
	int line_counter = 0;

	while(fgets(line, lineSize, file) != NULL){
		if(line_counter == line_number){
			line[strlen(line) - 1] = '\0';
			return line;
		} else {
			line_counter++;
		}
	}

	free(line);
	return NULL;
}

//Functie ce citeste dimensiunile imaginii, precum si 
//pixelii acesteia.
struct PGMImg read_image(const char *img_name, FILE* out){
	FILE* img_file = fopen(img_name, "r");
	if(img_file == NULL){
		printf("Eroare citire imagine!\n");
	}

	struct PGMImg img;
	int line_counter = 0;
	char line[100];
	char *token;
	const char s[] = " \n";

	//extragem header-ul fisierului .pgm
	char c;
	c = getc(img_file);
	int len = 0;
	while(c == 'P'){
		fputc(c, out);
		c = getc(img_file);
		while(c != '\n') {
			fputc(c, out);
			c = getc(img_file);
		}
		fputc('\n', out);
	}
	c = getc(img_file);

	while(c == '#'){
		fputc(c, out);
		c = getc(img_file);
		while(c != '\n') {
			fputc(c, out);
			c = getc(img_file);
		}
		fputc('\n', out);
	}

	//extragem dimensiunile matricii
	fscanf(img_file, "%d %d", &img.width, &img.height);
	fprintf(out, "%d ", img.width);
	fprintf(out, "%d\n", img.height);

	//extragem maximum_gray
	fscanf(img_file, "%d", &img.maximum_gray);
	fprintf(out, "%d\n", img.maximum_gray);

	//umplem matricea cu valorile pixelilor din fisier
	img.matrix = allocate_matrix(img.width + 2, img.height + 2);
	int i, j, aux;
	for(i = 0; i < img.height + 2; i++){
		for(j = 0; j < img.width + 2; j++){
			if(i == 0 || i == img.height + 1){
				img.matrix[i][j] = 0;
			} else if(j == 0 || j == img.width + 1){
				img.matrix[i][j] = 0;
			} else {
				fscanf(img_file, "%d", &aux);
				img.matrix[i][j] = aux;
			}
		}
	}

	return img;
}

//Functie ce extrage cei 8 vecini ai unui pixel.
void find_neighbours(int** matrix, int i, int j, int neighbours[8]){
	//NV
	neighbours[0] = matrix[i - 1][j - 1];
	//N
	neighbours[1] = matrix[i - 1][j];
	//NE
	neighbours[2] = matrix[i - 1][j + 1];
	//V
	neighbours[3] = matrix[i][j - 1];
	//E
	neighbours[4] = matrix[i][j + 1];
	//SV
	neighbours[5] = matrix[i + 1][j - 1];
	//S
	neighbours[6] = matrix[i + 1][j];
	//SE
	neighbours[7] = matrix[i + 1][j + 1];
}

//Functie ce aplica un filtru peste o matrice.
int** apply_filter(int **matrix, int width, int height, struct Filter filter){
	int** filtered_matrix = allocate_matrix(width + 2, height + 2);
	int neighbours[8];
	int i, j, m, n, new_val, counter, ok;
	int k;

	//parcurgem matricea fara borduri
	for(i = 1; i < height + 1; i++){
		for(j = 1; j < width + 1; j++){

			//identificam vecinii pixelului
			find_neighbours(matrix, i, j, neighbours);

			new_val = 0;
			counter = 0;
			ok = 0;

			//calculam noua valoare
			for(m = 0; m < 3; m++){
				for(n = 0; n < 3; n++){
					if(counter == 4 && ok!= 1){
						new_val += matrix[i][j] * filter.matrix[1][1];
						ok = 1;
					} else {
						new_val += neighbours[counter] * filter.matrix[m][n];
						counter++;
					}
				}
			}
			new_val /= filter.factor;
			new_val += filter.displacement;
				
			//daca valoarea pixelului depaseste intervalul [0...255]
			//ii dam cea mai apropiata valoare din interval
			if(new_val > 255){
				new_val = 255;
			} else if(new_val < 0){
				new_val = 0;
			}

			filtered_matrix[i][j] = new_val; 
		}
	}
	return filtered_matrix;
}

//Functie ce umple un vector cu o anumita valoare.
void set(int *a, int c, int size){
	int i;

	for (i = 0; i < size; ++ i){
		a[i] = c;
	}
}

//Functie ce intoarce numarul de vecini.
int neighbours_number(int *a, int size){
	int i, counter = 0;

	for(i = 0; i < size; i++){
		if(a[i] == 1){
			counter++;
		}
	}

	return counter;
}

int main(int argc, char *argv[]){
	char *topology;
	char *images;
	char *statistica;
	char *filter_name;
	int i, j;
	
	if(argc > 4){
		printf("Numar argumente prea mare!\n");
		return 0;
	} else {
		//primul paramteru este topologia
		topology = argv[1];
		//al doilea este fisierul cu iamgini
		images = argv[2];
		//al treilea este fisierul pentru statistica
		statistica = argv[3];
	}

	int rank;
	int nProcesses;
	int children_number = 0, portion = 0, counter = 0;

	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Request request;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	//citim numarul de imagini ce urmeaza sa fie procesate
	int images_nr = 0;
	FILE* images_in = fopen(images, "r");
	fscanf(images_in, "%d", &images_nr);

	int* neighbours = (int *)calloc(sizeof(int), nProcesses);
	int parent = -1;
	int limits[1];
	int null[0];
	set (null, 0,1);
	int* proc_lines = (int *)calloc(sizeof(int), nProcesses);
	int* recv_proc_lines = (int *)calloc(sizeof(int), nProcesses);
	int lines_per_proc = 0;

	int** new;
	int** to_be_filtered;
	int** filtered_matrix;
	int** final_matrix;

	//initializare vector vecini
	set_neighbours_array(neighbours, get_line(topology, rank), nProcesses);

	struct PGMImg img;
	struct Filter filter;

	int** temp_matrix;
	FILE* out;
	FILE* file_statistica = fopen(statistica, "w+");

	int curr_matrix_height;
	int tag, recv_tag;
	int k = 0;

	//Daca procesul are rank-ul 0, atunci el este radacina, adica initiatorul 
	//=> va trimite catre copiii toate infomatiile necesare prelucarii
	if(rank == 0){
		for(i = 0; i < images_nr; i++){
			char* line = get_line(images, i + 1);
			limits[0] = 1;
			limits[1] = 0;
			char *token;
			char *filter_name;
			char *image_in;
			char *image_out;
			int k = 0, width_sent = 0;
			counter = 0;
			int limits_receveid[1];

			//Extragem din fisierul de imagini numele filtrului,
			//imaginea de input si cea de output.
			const char s[] = " ";
			token = strtok(line, s);
			filter_name = token;
			while( token != NULL ) {
		      	token = strtok(NULL, s);
		      	if(k == 0){
		      		image_in = token;
		      	} else if(k == 1){
		      		image_out = token;
		      	}
		      	k++;
		   	}

		   	//Incepem sa trimitem imaginile.

		   	//Deschidem fisierul de out.
		   	out = fopen(image_out, "w+");

		   	//Citim imaginea ce urmeaza sa fie prelucrata.
		    img = read_image(image_in, out);

		    //Extragem numarul de copii al nodului.
			children_number = neighbours_number(neighbours, nProcesses);

			//Calculam dimensiunea pe care urmeaza sa o proceseze
			//fiecare copil.
			portion = img.height/children_number;
			limits[1] = portion;
			int p = 0;
			for(k = 0; k < nProcesses; k++){
				if(neighbours[k] == 1){
					
					//Pentru fiecare filtru avem un tag diferit.
					if(strcmp(filter_name, "sobel") == 0){
						tag = 1;
					} else if(strcmp(filter_name, "mean_removal") == 0){
						tag = 2;
					}

					//Trimitem copiilor latimea imaginii.
					MPI_Send(&img.width, 1, MPI_INT, k, tag, MPI_COMM_WORLD);

					//Trimitem copiilor limitele matricii.
					MPI_Send(limits, 2, MPI_INT, k, tag, MPI_COMM_WORLD);
					
					//Trimitem copiilor matricea linie cu linie.
					for(j = 0; j < img.height + 2; j++){
						if(j >= limits[0] - 1 && j <= limits[1] + 1){
							MPI_Send(img.matrix[j], img.width + 2, MPI_INT, k, tag, MPI_COMM_WORLD);
						}
					}
					limits[0] = limits[1] + 1;
					if((limits[1] + portion) < img.height && counter == (children_number - 2)){
						limits[1] = img.height;
					} else{
						limits[1] += portion;
					}
					counter++;
				}
			}

			//Alocam memorie pentru matricea in care vom
			//pastra imaginea filtrata.
			final_matrix = allocate_matrix(img.width  + 2, img.height);

			for(k = 0; k < nProcesses; k++){
				if(neighbours[k] == 1){
					//Asteptam sa primim limitele de la un copil.
					MPI_Recv(limits_receveid, 2, MPI_INT, k, tag, MPI_COMM_WORLD, &status);

					limits_receveid[0] -= 1;
					limits_receveid[1] -= 1;
					
					//Receptionam matricea prelucrata de copilul respectiv,
					//linie cu linie.
					for(j = 0; j < img.height; j++){
						if(j >= limits_receveid[0] && j <= limits_receveid[1]){
							MPI_Recv(final_matrix[j], img.width + 2, MPI_INT, k, tag, MPI_COMM_WORLD, &status);
						}
					}
				}
			}

			int m, n;
			//Scriem matricea in fisierul de output.
			for(m = 0; m < img.height; m++){
				for(n = 1; n <= img.width; n++){
					fprintf(out, "%d\n",final_matrix[m][n]);
				}
			}
		}
		
		//Trimitem tag-ul de terminare.
		for(k = 0; k < nProcesses; k++){
			if(neighbours[k] == 1){
				MPI_Send(null, 1, MPI_INT, k, 3, MPI_COMM_WORLD);
			}
		}

		//Primim vectorul cu liniile procesate.
		int* aux_recv_proc_lines =(int *)calloc(sizeof(int), nProcesses);
		for(k = 0; k < nProcesses; k++){
			if(neighbours[k] == 1){
				MPI_Recv(aux_recv_proc_lines, nProcesses, MPI_INT, k, 3, MPI_COMM_WORLD, &status);
				for(j = 0; j < nProcesses; j++){
					if(recv_proc_lines[j] == 0 && aux_recv_proc_lines[j] != 0){
						recv_proc_lines[j] = aux_recv_proc_lines[j];
					}
				}
			}
		}

		for(k = 0; k < nProcesses; k++){
			fprintf(file_statistica, "%d: %d\n", k, recv_proc_lines[k]);
		}

	} else if(rank != 0){

		lines_per_proc = 0;

		//Procesele care nu sunt initiator asteapta sa primeasca un mesaj de 
		//la alt proces care ii va deveni parinte.
		for(k = 0; k < images_nr; k++){
			int recv_limits[1], recv_limits_aux[1];
			int recv_width;
			struct Filter recv_filter;
			recv_filter.matrix = allocate_matrix(3, 3);

			//Primim latimea imaginii.
			MPI_Recv(&recv_width, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			recv_tag = status.MPI_TAG;
			parent = status.MPI_SOURCE;

			if(recv_tag == 1){
				recv_filter = readFilter("sobel");
			} else if(recv_tag == 2){
				recv_filter = readFilter("mean_removal");
			}

			//Primim limitele matricii.
			MPI_Recv(recv_limits, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			lines_per_proc += recv_limits[1] - recv_limits[0] + 1;

			//Alocam spatiu pentru portiunea din matricea initiala ce urmeaza 
			//sa fie primita si procesata.
			to_be_filtered = allocate_matrix(recv_width + 2,  recv_limits[1] - recv_limits[0] + 3);

			//Primim matricea linie cu linie
			for(i = 0; i <= recv_limits[1] - recv_limits[0] + 2; i++){
				MPI_Recv(to_be_filtered[i], recv_width + 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}

			limits[0] = 1;
			limits[1] = 0;

			//Daca nodul nu este frunza, dupa ce a receptionat mesajul, 
			//va imparti matricea in mod egal copiilor.
			counter = 0;
			int p = 0;
			if(neighbours_number(neighbours, nProcesses) != 1){
				//Calculam numaruld e copii.
				children_number = neighbours_number(neighbours,nProcesses) - 1;

				//Calculam bucata din matrice ce urmeaza sa fie procesata 
				// de copii.
				portion = (recv_limits[1] - recv_limits[0] +1)/children_number;

				limits[1] = portion;
				for(i = 0; i < nProcesses; i++){
					if(neighbours[i] == 1 && i != parent){
						
						//Trimitem copiilor latimea matricii.
						MPI_Send(&recv_width, 1, MPI_INT, i, recv_tag, MPI_COMM_WORLD);

						//Trimitem copiilor limitele matricii.
						MPI_Send(limits, 2, MPI_INT, i, recv_tag, MPI_COMM_WORLD);

						//Trimitem matricea linie cu linie copiilor/
						for(j = 0; j <= recv_limits[1] - recv_limits[0] + 2; j++){
							if(j >= limits[0] - 1 && j <= limits[1] + 1){
								MPI_Send(to_be_filtered[j], recv_width + 2, MPI_INT, i, recv_tag, MPI_COMM_WORLD);

							}
						}

						limits[0] = limits[1] + 1;
						if((limits[1] + portion) < (recv_limits[1] - recv_limits[0] + 2) && counter == (children_number - 2)){
							limits[1] = recv_limits[1] - recv_limits[0] + 1;
						} else{
							limits[1] += portion;
						}
						counter++;					
					}
				}

				//Alocam memorie pentru matricea in care vom pastra portiunea 
				//din matricea initiala filtrata. 
				filtered_matrix = allocate_matrix(recv_width  + 2,  recv_limits[1] - recv_limits[0] + 1);

				int p = 0;
				for(i = 0; i <  nProcesses; i++){
					if(neighbours[i] == 1 && i != parent){
						//Primeste limitele de la copii.
						MPI_Recv(recv_limits_aux, 2, MPI_INT, i, recv_tag, MPI_COMM_WORLD, &status);
						recv_limits_aux[0] -= 1;
						recv_limits_aux[1] -= 1;

						for(j = 0; j <= recv_limits[1] - recv_limits[0]; j++){
							if(j >= recv_limits_aux[0] && j <= recv_limits_aux[1]){
								MPI_Recv(filtered_matrix[j], recv_width + 2, MPI_INT, i, recv_tag, MPI_COMM_WORLD, &status);
							}
						}
					}
				}

				//Trimitem parintelui limitele matricii.
				MPI_Send(recv_limits, 2, MPI_INT, parent, recv_tag, MPI_COMM_WORLD);

				//Trimitem parintelui matricea.
				for(i = 0; i <  nProcesses; i++){
					if(neighbours[i] == 1 && i == parent){
						for(j = 0; j <= recv_limits[1] - recv_limits[0]; j++){
							MPI_Send(filtered_matrix[j], recv_width + 2, MPI_INT, i, recv_tag, MPI_COMM_WORLD);
						}
					}
				}

			} else if(neighbours_number(neighbours, nProcesses) == 1){
					//Daca nodul este frunza, acesta trebuie sa aplice filtrul //pe portiunea de matrice primita.
					
					//Alocam spatiu pentru ,atricea ce va fi filtrata.
					new = allocate_matrix(recv_width + 2,  recv_limits[1] - recv_limits[0] + 3);
					new = apply_filter(to_be_filtered, recv_width, recv_limits[1] - recv_limits[0] + 1, recv_filter);

					//Trimitem parintelui limitele matricii.
					MPI_Send(recv_limits, 2, MPI_INT, parent, recv_tag, MPI_COMM_WORLD);

					int p = 0;
					//Trimitem matricea linie cu linie parintelui.
					for(j = 1; j < recv_limits[1] - recv_limits[0] + 2; j++){
						MPI_Send(new[j], recv_width + 2, MPI_INT, parent, recv_tag, MPI_COMM_WORLD);
					}
				}
			}

			//Receptionam tag-ul de terminare.
			int recv_null[0];
			MPI_Recv(recv_null, 1, MPI_INT, parent, 3, MPI_COMM_WORLD, &status);

			//Daca procesul este intermediar se da mai departe copiilor.
			if(neighbours_number(neighbours, nProcesses) != 1){
				for(i = 0; i < nProcesses; i++){
					if(neighbours[i] == 1 && i != parent){
						MPI_Send(null, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
					}
				}
			} else if(neighbours_number(neighbours, nProcesses) == 1){
				//Daca procesul este frunza trimite parintelui numarul de linii procesate.
				for(i = 0; i < nProcesses; i++){
					if(i == rank){
						proc_lines[i] = lines_per_proc;
					}
				}
				MPI_Send(proc_lines, nProcesses, MPI_INT, parent, 3, MPI_COMM_WORLD);
			}

			//Receptionam liniile procesate si le trimitem mai departe 
			//parintelui.
			int* aux_proc_lines = (int *)calloc(sizeof(int), nProcesses);
			if(neighbours_number(neighbours, nProcesses) != 1){
				for(i = 0; i < nProcesses; i++){
					if(neighbours[i] == 1 && i != parent){
						MPI_Recv(aux_proc_lines, nProcesses, MPI_INT, i, 3, MPI_COMM_WORLD, &status);
						for(j = 0; j < nProcesses; j++){
							if(proc_lines[j] == 0 && aux_proc_lines[j] != 0){
								proc_lines[j] = aux_proc_lines[j];
							}
						}
					}
				}
				MPI_Send(proc_lines, nProcesses, MPI_INT, parent, 3, MPI_COMM_WORLD);
			}
		}

	MPI_Finalize();
	return 0;
}
