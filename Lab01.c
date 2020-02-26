#include <stdio.h>
#include <stdlib.h>

//--------------------------------------------------

const int FIELD_FOR_COMPLEX = 10;

typedef struct ComplexNumber {
    float re;
    float im;
} comp;

void CreateComp(comp *tmp, float re, float im) {
    tmp->re = re;
    tmp->im = im;
}

void FreeComp(comp *tmp) {
    free(tmp);
}

void CopyComp(comp *from, comp *to) {
    to->re = from->re;
    to->im = from->im;
}

void CompSumm(comp *z1, comp *z2, comp *res) {
    float a = z1->re;
    float b = z1->im;
    float c = z2->re;
    float d = z2->im;
    res->re = a + c;
    res->im = b + d;
}

void CompDist(comp *z1, comp *z2, comp *res) {
    float a = z1->re;
    float b = z1->im;
    float c = z2->re;
    float d = z2->im;
    res->re = a - c;
    res->im = b - d;
}

void CompMult(comp *z1, comp *z2, comp *res) {
    float a = z1->re;
    float b = z1->im;
    float c = z2->re;
    float d = z2->im;
    res->re = (a * c) - (b * d);
    res->im = (a * d) + (b * c);
}

void CompDiv(comp *z1, comp *z2, comp *res) {
    float a = z1->re;
    float b = z1->im;
    float c = z2->re;
    float d = z2->im;
    res->re = ((a * c) + (b * d)) / (c * c + d * d);
    res->im = (b * c) - (a * d) / (c * c + d * d);
}

void IntoComp(comp *tmp, float num) {
    tmp->re = num;
    tmp->im = 0;
}

void PrintComp(comp *tmp) {
    if (tmp->re == 0) {
        if (tmp->im == 0) {
            printf("%10d", 0);
        } else {
            printf("%10.1fi", tmp->im);
        }
    } else {
        if (tmp->im == 0) {
            printf("%10.1f", tmp->re);
        } else if (tmp->im > 0) {
            printf("%4.1f%s%4.1fi", tmp->re, "+", tmp->im);
        } else {
            printf("%4.1f%4.1fi", tmp->re, tmp->im);
        }
    }
}

//--------------------------------------------------
typedef struct Matrix {
    int rows;
    int cols;
    comp **content;
} matr;

void FreeMatrix(matr *tmp) {
    for (int row = 0; row < tmp->rows; ++row) {
        free(tmp->content[row]);
    }
    tmp->content = NULL;
}

void CreateMatr(matr *tmp, int rows, int cols, comp *content) {
    if (tmp != NULL)
        FreeMatrix(tmp);
    tmp->rows = rows;
    tmp->cols = cols;
    if (content != NULL) {
        int cnt = 0;
        tmp->content = (comp **) malloc(sizeof(comp *) * rows);
        for (int row = 0; row < rows; ++row) {
            tmp->content[row] = (comp *) malloc(sizeof(comp) * cols);
            for (int col = 0; col < cols; ++col) {
                CopyComp(content + cnt, (tmp->content[row] + col));
                ++cnt;
            }
        }
    } else {
        comp *zero = (comp *) malloc(sizeof(comp));
        IntoComp(zero, 0);
        tmp->content = (comp **) malloc(sizeof(comp *) * rows);
        for (int row = 0; row < rows; ++row) {
            tmp->content[row] = (comp *) malloc(sizeof(comp) * cols);
            for (int col = 0; col < cols; ++col) {
                CopyComp(zero, (tmp->content[row] + col));
            }
        }
        FreeComp(zero);
    }
}

void MatrSumm(matr *a1, matr *a2, matr *res) {
    if ((a1->rows == a2->rows) && (a1->cols == a2->cols)) {
        comp *content = (comp *) malloc(sizeof(comp) * a1->rows * a1->cols);
        int cnt = 0;
        for (int row = 0; row < a1->rows; ++row) {
            for (int col = 0; col < a1->cols; ++col) {
                CompSumm(a1->content[row] + col, a2->content[row] + col, content + cnt);
                ++cnt;
            }
        }
        CreateMatr(res, a1->rows, a1->cols, content);
        FreeComp(content);
    } else {
        printf("\nThis matrixes cannot be summarized\n");
    }
}

void MatrMult(matr *a1, matr *a2, matr *res) {
    if (a1->cols == a2->rows) {
        comp *tmp = (comp *) malloc(sizeof(comp));
        comp *final = (comp *) malloc(sizeof(comp));
        IntoComp(tmp, 0);
        IntoComp(final, 0);
        comp *content = (comp *) malloc(sizeof(comp) * a1->rows * a2->cols);
        int cnt = 0;
        for (int row = 0; row < a1->rows; ++row) {
            for (int col = 0; col < a2->cols; ++col) {
                for (int flag = 0; flag < a1->cols; ++flag) {
                    CompMult(a1->content[row] + flag, a2->content[flag] + col, tmp);
                    CompSumm(tmp, final, final);
                }
                CopyComp(final, (content + cnt));
                ++cnt;
                IntoComp(tmp, 0);
                IntoComp(final, 0);
            }
        }
        CreateMatr(res, a1->rows, a2->cols, content);
        FreeComp(tmp);
        FreeComp(final);
        FreeComp(content);
    } else {
        printf("\nThis matrixes cannot be multiplied\n");
    }
}

void PrintMatr(matr *tmp) {

    int rows = tmp->rows;
    int cols = tmp->cols;
    printf("\nSize%d%s%d%s", rows, "x", cols, "\n ");
    for (int cnt = 0; cnt < cols * FIELD_FOR_COMPLEX + cols * 3; ++cnt) {
        printf("-");
    }
    printf("\n");
    for (int row = 0; row < rows; ++row) {
        printf("%s", "| ");
        for (int col = 0; col < cols; ++col) {
            PrintComp(tmp->content[row] + col);
            printf(" | ");
        }
        printf("\n ");
        for (int cnt = 0; cnt < cols * FIELD_FOR_COMPLEX + cols * 3; ++cnt) {
            printf("-");
        }
        printf("\n");
    }
}

void MatrTranspose(matr *a, matr *res) {
    CreateMatr(res, a->cols, a->rows, NULL);
    for (int row = 0; row < a->rows; ++row) {
        for (int col = 0; col < a->cols; ++col) {
            CopyComp(*(a->content + row) + col, *(res->content + col) + row);
        }
    }
}

void Interface() {
    const int MAX_STR_LENGTH = 20;
    int Arows = 2;
    int Acols = 2;
    int Brows = 2;
    int Bcols = 2;
    comp *Acontent = malloc(sizeof(comp) * 4);
    comp *Bcontent = malloc(sizeof(comp) * 4);
    comp *tmp = malloc(sizeof(comp));
    for (int i = 0; i < 4; ++i) {
        CreateComp(tmp, i, 0);
        CopyComp(tmp, (Acontent + i));
        CopyComp(tmp, (Bcontent + i));
        PrintComp(Acontent + i);
        PrintComp(Bcontent + i);
        printf("\n");
    }
    matr *A = malloc(sizeof(matr));
    CreateMatr(A, Arows, Acols, Acontent);
    matr *B = malloc(sizeof(matr));
    CreateMatr(B, Brows, Bcols, Bcontent);
    printf("%s", "Default values:\n");
    printf("%s", "A:");
    PrintMatr(A);
    printf("%s", "B:");
    PrintMatr(B);
    int inp = 9;
    while (1) {
        printf("%s", "Available commands:"
                     "\n1 - Change matrix A;"
                     "\n2 - Change matrix B;"
                     "\n3 - Show matrix A"
                     "\n4 - Show matrix B"
                     "\n5 - Show A+B"
                     "\n6 - Show A*B"
                     "\n7 - Show transposed A"
                     "\n8 - Show transposed B"
                     "\n9 - Exit\n");
        int inp;
        scanf("%d", &inp);
        if (inp == 1) {
            int rows = 0;
            int cols = 0;
            printf("How many rows and collums should matrix A have?\nFormat: rows collums\n>");
            scanf("%d %d", &rows, &cols);
            while (rows < 1 || cols < 1) {
                printf("Rows and collums must be natural number!\n");
                scanf("%d %d", &rows, &cols);
            }
            comp *content = (comp *) malloc(sizeof(comp) * rows * cols);
            comp *tmp = (comp *) malloc(sizeof(comp));
            printf("Input numbers that you want to be in matrix A:\n");
            for (int i = 0; i < rows * cols; ++i) {
                printf("Real and Imaginary parts of complex number %d%s", i, " = ");
                float tmpRe = 0;
                float tmpIm = 0;
                scanf("%f %f", &tmpRe, &tmpIm);
                CreateComp(tmp, tmpRe, tmpIm);
                CopyComp(tmp, (content + i));
            }
            CreateMatr(A, rows, cols, content);
            FreeComp(tmp);
            FreeComp(content);
        } else if (inp == 2) {
            int rows = 0;
            int cols = 0;
            printf("How many rows and collums should matrix B have?\nFormat: rows collums\n>");
            scanf("%d %d", &rows, &cols);
            while (rows < 1 || cols < 1) {
                printf("Rows and collums must be natural number!\n");
                scanf("%d %d", &rows, &cols);
            }
            comp *content = (comp *) malloc(sizeof(comp) * rows * cols);
            comp *tmp = (comp *) malloc(sizeof(comp));
            printf("Input numbers that you want to be in matrix B:\n");
            for (int i = 0; i < rows * cols; ++i) {
                printf("Real and Imaginary parts of complex number %d%s", i, " = ");
                float tmpRe = 0;
                float tmpIm = 0;
                scanf("%f %f", &tmpRe, &tmpIm);
                CreateComp(tmp, tmpRe, tmpIm);
                CopyComp(tmp, (content + i));
            }
            CreateMatr(B, rows, cols, content);
            FreeComp(tmp);
            FreeComp(content);
        } else if (inp == 3) {
            printf("Matrix A\n");
            PrintMatr(A);
        } else if (inp == 4) {
            printf("Matrix B\n");
            PrintMatr(B);
        } else if (inp == 5) {
            matr *res = (matr *) malloc(sizeof(matr));
            MatrSumm(A, B, res);
            printf("A+B\n");
            PrintMatr(res);
            FreeMatrix(res);
        } else if (inp == 6) {
            matr *res = (matr *) malloc(sizeof(matr));
            MatrMult(A, B, res);
            printf("A*B\n");
            PrintMatr(res);
            FreeMatrix(res);
        } else if (inp == 7) {
            matr *tmp = malloc(sizeof(matr));
            MatrTranspose(A, tmp);
            PrintMatr(tmp);
            FreeMatrix(tmp);
        } else if (inp == 8) {
            matr *tmp = malloc(sizeof(matr));
            MatrTranspose(A, tmp);
            PrintMatr(tmp);
            FreeMatrix(tmp);
        } else if (inp == 9) {
            printf("\n\nShutting down the interface.\n\n");
            free(Acontent);
            free(Bcontent);
            FreeComp(tmp);
            FreeMatrix(A);
            FreeMatrix(B);
            break;
        } else {
            printf("Unknown command\n\n");
        }
    }
}

void ShowInfo(matr *A, matr *B, matr *res) {
    printf("%s", "\nMatrix A:");
    PrintMatr(A);
    printf("%s", "\nMatrix B:");
    PrintMatr(B);

    printf("%s", "\nA+B");
    MatrSumm(A, B, res);
    if (res->rows > 0 && res->cols > 0) {
        PrintMatr(res);
    }
    printf("%s", "\nA*B");
    MatrMult(A, B, res);
    if (res->rows > 0 && res->cols > 0) {
        PrintMatr(res);
    }
    printf("%s", "\nTransposed A");
    MatrTranspose(A, res);
    PrintMatr(res);
    printf("%s", "\nTransposed B");
    MatrTranspose(B, res);
    PrintMatr(res);
    FreeMatrix(res);
}

void Tests() {
    matr *A = (matr *) malloc(sizeof(matr));
    matr *B = (matr *) malloc(sizeof(matr));
    comp *tmp = (comp *) malloc(sizeof(comp));

    comp *fcont = (comp *) malloc(sizeof(comp) * 9);
    for (int i = 0; i < 9; ++i) {
        CopyComp(tmp, (fcont + i));
    }
    int i = 0;
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            if (row == col) {
                IntoComp(tmp, 1);
            } else {
                IntoComp(tmp, 0);
            }
            CopyComp(tmp, (fcont + i));
            ++i;
        }
    }
    CreateMatr(A, 3, 3, fcont);

    fcont = (comp *) malloc(sizeof(comp) * 12);
    for (int i = 0; i < 12; ++i) {
        CreateComp(tmp, i, i);
        CopyComp(tmp, (fcont + i));
    }
    CreateMatr(B, 3, 4, fcont);
    FreeComp(fcont);

    matr *res = (matr *) malloc(sizeof(matr));
    res->rows = 0;
    res->cols = 0;
    printf("%s", "\nTEST #1\n");
    ShowInfo(A,B,res);

    fcont = (comp *) malloc(sizeof(comp) * 9);
    for (int i = 0; i < 9; ++i) {
        IntoComp(tmp, i);
        CopyComp(tmp, fcont + i);
    }
    CreateMatr(A, 3, 3, fcont);
    CreateMatr(B, 3, 3, fcont);

    res = (matr *) malloc(sizeof(matr));
    res->rows = 0;
    res->cols = 0;
    printf("%s", "\nTEST #2\n");
    ShowInfo(A,B,res);

    FreeComp(fcont);
    fcont = (comp *) malloc(sizeof(comp) * 3);
    for (int i = 0; i < 3; ++i) {
        CreateComp(tmp, i, i);
        CopyComp(tmp, fcont + i);
    }
    CreateMatr(A, 3, 1, fcont);
    FreeComp(fcont);
    fcont = (comp *) malloc(sizeof(comp));
    IntoComp(tmp, 1);
    CopyComp(tmp, fcont);
    CreateMatr(B, 1, 1, fcont);

    res = (matr *) malloc(sizeof(matr));
    res->rows = 0;
    res->cols = 0;
    printf("%s", "\nTEST #3\n");
    ShowInfo(A,B,res);

    FreeComp(fcont);
    fcont = (comp *) malloc(sizeof(comp) * 4);
    for (int i = 0; i < 4; ++i) {
        IntoComp(tmp, i);
        CopyComp(tmp, fcont + i);
    }
    CreateMatr(A, 2, 2, fcont);
    FreeComp(fcont);
    fcont = (comp *) malloc(sizeof(comp) * 9);
    for (int i = 0; i < 9; ++i) {
        IntoComp(tmp, i);
        CopyComp(tmp, fcont + i);
    }
    CreateMatr(B, 3, 3, fcont);

    res = (matr *) malloc(sizeof(matr));
    res->rows = 0;
    res->cols = 0;
    printf("%s", "\nTEST #4\n");
    ShowInfo(A,B,res);

    FreeComp(fcont);
    FreeComp(tmp);
    FreeMatrix(A);
    FreeMatrix(B);
}

//--------------------------------------------------
int main() {
    char symb;
    while (1) {
        printf("%s", "\nLab01.\nTo run tests enter 't',\nFor interface enter 'i'\n to exit enter 'x'.\n");
        scanf("%c", &symb);
        if (symb == 'i') {
            Interface();
        } else if (symb == 't') {
            Tests();
        } else if (symb == 'x') {
            printf("\n\nShutting down the program.\n\n");
            break;
        }
        scanf("%c", &symb);
    }
    return 0;
}