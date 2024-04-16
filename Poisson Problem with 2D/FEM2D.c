#include <stdio.h>
#include <stdlib.h>

void read_mesh();
int main(int argc, char *argv[])
{
  
  return 0;
}

void read_mesh() {
  FILE* mesh;
  mesh = fopen("mesh.dat","r");
  if (mesh == NULL) {
    perror("Error while opening the file.\n");
    exit(EXIT_FAILURE);
  }
  while (fscanf(ifp, "%d %d", username, &score) != EOF) {
    fprintf(ofp, "%s %d\n", username, score+10);
  }
  fclose(mesh);
}
