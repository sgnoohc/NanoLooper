#include<stdio.h>
#include<dirent.h>


int main(argc, argv)
int argc;
char *argv[];
{
    DIR *d;
    char *filename = "commands.txt";
    FILE *output = fopen(filename, "w");
    struct dirent *dir;
    d = opendir("/home/users/joytzphysics/ttbar/.");
    int n = 0;
    if (output == NULL)
    {
        printf("Error opening the file %s", filename);
        return -1;
    }
    if (d)
    {
        while ((dir = readdir(d)) != NULL)
        {
            fprintf(output,"./runNanoLooper --input /home/users/joytzphysics/ttbar/%s --output /home/users/joytzphysics/NanoLooper/ttbaroutput/output_%d.root --scale1fb 0.066\n", dir->d_name, n);
            n++;
        }
        closedir(d);
    }
    fclose(output);

    return 0;
}

    
