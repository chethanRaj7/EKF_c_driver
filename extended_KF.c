#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float voltage[100], current[100], Temp[100], curr_load[100], volt_load[100], time[100], cap[100];

void choose_variable(int variable, float value, int index){
    switch (variable) {
    case 1:
        voltage[index] = value;
        printf("Voltage: %f\n", value);
        break;
    case 2:
        current[index] = value;
        printf("current: %f\n", value);
        break;
    case 3:
        Temp[index] = value;
        printf("temp: %f\n", value);
        break;
    case 4:
        curr_load[index] = value;
        printf("curr_load: %f\n", value);
        break;
    case 5:
        volt_load[index] = value;
        printf("vol_load: %f\n", value);
        break;
    case 6:
        time[index] = value;
        printf("time: %f\n", value);
        break;
    case 7:
        cap[index] = value;
        printf("cap: %f\n", value);
        break;
    default:
        printf("default value");
}
}

int count_lines(FILE* file)
{
    char buf[65536];
    int counter = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, 65536, file);
        if (ferror(file))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    return counter;
}

void importdata(FILE *fp, int num_lines)
{
    char chunk[124];
    char letter[20];
    int num_variables = 8;
    int temp_index = 0;

    //int num_lines = count_lines(fp);
    //printf("Num of lines: %d\n", num_lines);
    fseek(fp, 0, SEEK_SET);

    while(temp_index< num_lines){
        fgets(chunk, sizeof(chunk), fp) != NULL;
        //printf("%s\n", chunk);
        
        //populate each variable
        int j = 1;
        int temp_i = 0;
        //Get each variable.
        for(int i = 0; i< sizeof(chunk); i++){

            if (chunk[i] == '\0'){
                char *endptr;                                        // Pointer to track invalid characters
                float temp_value = strtof(letter, &endptr);          // Convert string to float

                //float temp_value = strtod(temp_value, NULL);
                printf("populating: %d\n", j);
                choose_variable(j, temp_value, temp_index);           //populate the values to index.
                temp_i = 0;                                           //reset the index after update

                break;
            }
            else if (chunk[i] != ','){
                letter[temp_i] = chunk[i];
                temp_i++;
            }
            else if(chunk[i] == ','){
                char *endptr;                                        // Pointer to track invalid characters
                float temp_value = strtof(letter, &endptr);          // Convert string to float

                printf("populating: %d\n", j);
                choose_variable(j, temp_value, temp_index);           //populate the values to index.
                j++;
                temp_i = 0;                                           //reset the index after update
                continue;
            }
        }
        if (temp_index == 10){
            break;
        }
        temp_index++;
    }
    fclose(fp);
}


void main()
{
    float xHat[] = {1.0, 0.0, 0.0};
    float sigmaX[3][3] = {{1e-5, 0.0, 0.0},
                        {0.0, 1e-5, 0.0},
                        {0.0, 0.0, 1e-5}
                        };

    float socVar = 1e1;
    float iRC1Var = 1e1;
    float iRC2Var = 1e1;

    float sigmaW[3][3] = {{socVar, 0.0, 0.0},
                        {0.0, iRC1Var, 0.0},
                        {0.0, 0.0, iRC2Var}
                        };

    float sigmaV = 1e-3;
    float socStore[490] = {0};


    FILE *fp = fopen("00001.csv", "r");
    if(fp == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    int num_lines = count_lines(fp);
    

    // float *voltage = (float *)malloc(num_lines * sizeof(float));
    // float *current = (float *)malloc(num_lines * sizeof(float));
    // float *Temp = (float *)malloc(num_lines * sizeof(float));
    // float *curr_load = (float *)malloc(num_lines * sizeof(float));
    // float *volt_load = (float *)malloc(num_lines * sizeof(float));
    // float *time = (float *)malloc(num_lines * sizeof(float));
    // float *cap = (float *)malloc(num_lines * sizeof(float));

    importdata(fp, num_lines);

}