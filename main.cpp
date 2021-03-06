//zapisuje do pliku analityczna i kmb



#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include "calerf.h"

using namespace std;


//stale z polecenia
const double D = 1.0;
const double lambdaKmb = 0.4;
const double R = 1.0;
const double A = 10.0;
const double tMax = 2.0;
fstream daneAnalityczny;

//ilosc poziomow czasowych
int sizeKmb, xSize;

void analityczne_KMB(double **W, double dt, double h)
{
    double tt = dt;
    double dh = 0;

//dla jednego poziomu czasowego
    for(int i = 1; i < sizeKmb; i++)
    {
        dh = 1;
//dla kazdego x
        for(int j = 0; j < xSize; j++)
        {
            W[i][j] = 1.0-(((double)R/(double)dh)*erfc((double)(dh-R)/(2.0*sqrt(D*tt))));
            dh += h;
        }
        tt += dt;
    }
}

//warunek brzegowy. 
double war_prawy(double t) {
    return 1.0 - (R / (R + A)) * (1 - erfl(A / (2.0 * sqrt(D * t))));
}

//funkcja KMB
void fun_KMB(double **W, double h) {
    double dh;

    for (int i = 1; i < sizeKmb; i++) {
        dh = R + h;
//wzor z wykladu. 
//mowil ze trzeab sprawdzac od lambd mniejszych az do takich najblizszej lambdzie
//a im h sie zmniejsza tym lambda bedzie blizsza 1. 
//dla srodkowego nie robimy tak bo nie ma przejscia miedzy wezlami
        for (int j = 1; j < xSize - 1; j++) {
            W[i][j] = W[i - 1][j - 1] * lambdaKmb * (1.0 - h / dh)
                      + W[i - 1][j] * (1.0 - (2.0 * lambdaKmb))
                      + W[i - 1][j + 1] * lambdaKmb * (1.0 + h / dh);
            dh += h;
        }
    }
}

int main()
{
    double bladCrankNicolson, bladKmb, maxBladCrankNicolson, maxBladKmb, h = 0.5;

    for (; h > 0.01; h -= 0.01) {
        //liczymy ile poziomow czasowych
        xSize = A / h;

        double **Analityczne_K, **KMB;
        double dt_kmb = (lambdaKmb * h * h) / D;
        //krok czasowy z wzory wziety ( lambda = h*h * dt)/D
        //to jest ten krok dt dla kazdego 

        sizeKmb = tMax / dt_kmb;

        Analityczne_K = new double *[sizeKmb];
        KMB = new double *[sizeKmb];

        for (int i = 0; i < sizeKmb; i++) {
            Analityczne_K[i] = new double[xSize];
            KMB[i] = new double[xSize];
        }

        //warunek poczatkowy
        for (int i = 0; i < xSize; i++) {
            KMB[0][i] = 1.0;
        }

         //warunek brzegowy 1
        for (int i = 1; i < sizeKmb; i++) {
            KMB[i][0] = 0.0;
        }

        double tt = dt_kmb;

        for (int i = 1; i < sizeKmb; i++) {
            KMB[i][xSize - 1] = war_prawy(tt);
            tt += dt_kmb;
        }

        analityczne_KMB(Analityczne_K, dt_kmb, h);
        fun_KMB(KMB, h);

        if (h < 0.1001 && h > 0.0999) {
            int index = 150; //50 lub 100 lub 150
            int index2 = index * 2 + (index/50)*25;

            double dh = R + h;
            daneAnalityczny.open("daneAnalityczny.txt", fstream::out);

            for (int i = 0; i < xSize; i++) {
                daneAnalityczny << dh << " " << " " <<
                                Analityczne_K[index2][i] << " " << KMB[index2][i] <<  endl;
                dh += h;
                cout.width(4);
                cout << dh << " ";
                cout.width(15);
                cout << setprecision(10) << Analityczne_K[index2][i] << " ";
                cout.width(15);
                cout << setprecision(10) << KMB[index2][i] << " " << endl;
            }
            daneAnalityczny.close();
        }
    }
    return 0;
}








/////wsyzstko dla kmb



#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include "calerf.h"

using namespace std;


//stale z polecenia
const double D = 1.0;
const double lambdaKmb = 0.4;
const double R = 1.0;
const double A = 10.0;
const double tMax = 2.0;
fstream daneAnalityczny, daneErrKmb, daneDlaTmax;

//ilosc poziomow czasowych
int sizeKmb, xSize;

void analityczne_KMB(double **W, double dt, double h)
{
    double tt = dt;
    double dh = 0;

//dla jednego poziomu czasowego
    for(int i = 1; i < sizeKmb; i++)
    {
        dh = 1;
//dla kazdego x
        for(int j = 0; j < xSize; j++)
        {
            W[i][j] = 1.0-(((double)R/(double)dh)*erfc((double)(dh-R)/(2.0*sqrt(D*tt))));
            dh += h;
        }
        tt += dt;
    }
}

//warunek brzegowy. 
double war_prawy(double t) {
    return 1.0 - (R / (R + A)) * (1 - erfl(A / (2.0 * sqrt(D * t))));
}

//funkcja KMB
void fun_KMB(double **W, double h) {
    double dh;

    for (int i = 1; i < sizeKmb; i++) {
        dh = R + h;
//wzor z wykladu. 
//mowil ze trzeab sprawdzac od lambd mniejszych az do takich najblizszej lambdzie
//a im h sie zmniejsza tym lambda bedzie blizsza 1. 
//dla srodkowego nie robimy tak bo nie ma przejscia miedzy wezlami
        for (int j = 1; j < xSize - 1; j++) {
            W[i][j] = W[i - 1][j - 1] * lambdaKmb * (1.0 - h / dh)
                      + W[i - 1][j] * (1.0 - (2.0 * lambdaKmb))
                      + W[i - 1][j + 1] * lambdaKmb * (1.0 + h / dh);
            dh += h;
        }
    }
}

int main()
{
    double bladCrankNicolson, bladKmb, maxBladCrankNicolson, maxBladKmb, h = 0.5;
    daneDlaTmax.open("daneDlaTmax.txt", fstream::out);
    for (; h > 0.01; h -= 0.01) {
        //liczymy ile poziomow czasowych
        xSize = A / h;

        double **Analityczne_K, **KMB;
        double dt_kmb = (lambdaKmb * h * h) / D;
        //krok czasowy z wzory wziety ( lambda = h*h * dt)/D
        //to jest ten krok dt dla kazdego 

        
        sizeKmb = tMax / dt_kmb; //ile iteracji (Tmax / dtkmb)

        //Tworzymy macierze dla KMB, jedne to po prostu analityczne
        //druga liczona metoda KMB
        Analityczne_K = new double *[sizeKmb];
        KMB = new double *[sizeKmb];

        for (int i = 0; i < sizeKmb; i++) {
            Analityczne_K[i] = new double[xSize];
            KMB[i] = new double[xSize];
        }

        //warunek poczatkowy
        //z polecenia
        for (int i = 0; i < xSize; i++) {
            KMB[0][i] = 1.0;
        }

         //warunek brzegowy 1
         //z polecenia
        for (int i = 1; i < sizeKmb; i++) {
            KMB[i][0] = 0.0;
        }

        //pomocnicze
        double tt = dt_kmb;

        //zastapienie drugiego warunku tym:
        for (int i = 1; i < sizeKmb; i++) {
            KMB[i][xSize - 1] = war_prawy(tt);
            tt += dt_kmb;
        }

        //rozwiazanie analityczne
        analityczne_KMB(Analityczne_K, dt_kmb, h);
        //policzenie za pomoca kmb
        fun_KMB(KMB, h);

        //dla h = 0.1
        if (h < 0.1001 && h > 0.0999) {
            int index = 50; //50 lub 100 lub 150 w zaleznosci od 0.5. 1. 1,5
            int index2 = index * 2 + (index/50)*25;

            //o krok
            //-----tutaj wykres od xow, od 1 do 10. wartosci, na jednym poziomie czasowym
            double dh = R + h;
            daneAnalityczny.open("daneAnalityczny.txt", fstream::out);

            for (int i = 0; i < xSize; i++) {
                daneAnalityczny << dh << " " << " " <<
                                Analityczne_K[index2][i] << " " << KMB[index2][i] <<  endl;
                dh += h;
                cout.width(4);
                cout << dh << " ";
                cout.width(15);
                cout << setprecision(10) << Analityczne_K[index2][i] << " ";
                cout.width(15);
                cout << setprecision(10) << KMB[index2][i] << " " << endl;
            }
            daneAnalityczny.close();

            double errCrankN, errKmb, maxErrCrankN, maxErrKmb;

            //dla kazdego x, max blad w czasie. 
            daneErrKmb.open("daneErrKmb.txt", fstream::out);
            for (int i = 1; i < sizeKmb; i++) {
                maxErrKmb = 0.0;
                for (int j = 0; j < xSize; j++) {
                    errKmb = fabs(Analityczne_K[i][j] - KMB[i][j]);
                    if (errKmb > maxErrKmb)
                        maxErrKmb = errKmb;
                }
                daneErrKmb << i * dt_kmb << " " << maxErrKmb << endl;
            }

            daneErrKmb.close();
        }
        maxBladCrankNicolson = 0.0;
        maxBladKmb = 0.0;
        //Maksymalna wartosc bezwzglednego bledu obserwacji dla tMax w funkcji kroku przestrzennnego h
        for (int i = 0; i < xSize; i++) {
            // bladCrankNicolson = fabs(AnalCrankNicolson[sizeCrankNicolson - 1][i] - CrankNicolson[sizeCrankNicolson - 1][i]);
            // if (bladCrankNicolson > maxBladCrankNicolson)
            //     maxBladCrankNicolson = bladCrankNicolson;

            bladKmb = fabs(Analityczne_K[sizeKmb - 1][i] - KMB[sizeKmb - 1][i]);
            if (bladKmb > maxBladKmb)
                maxBladKmb = bladKmb;

        }
        //daneDlaTmax << log10(h) << " " << log10(maxBladCrankNicolson) << " " << log10(maxBladKmb) << endl;
        daneDlaTmax << log10(h) << " " << log10(maxBladKmb) << endl;
        
    }
    daneDlaTmax.close();
    return 0;
}
