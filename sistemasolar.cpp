//---------------------------------------------------------------------
//                                                                    |
//       ALGORITMO DE VERLET PARA EL SISTEMA SOLAR (9 PLANETAS)       |
//                                                                    |
//                            Objetivos:                              |
//                                                                    |
//  1.- Hacer vídeo de las órbitas                                    |
//  2.- Plot energía E(t)                                             |
//  3.- Calcular el período de rotación de los planetas               |
//  4.- Comprobar que el sistema es estable frente a perturbaciones.  |
//                                                                    |
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//                                                                    |
//                       CONDICIONES INICIALES                        |
//   windows2universe.org/our_solar_system/planets_orbits_table.html  |
//                                                                    |
//    i     Planeta      T(años)    Vel. orbital (m/s)    Masa (kg)   |
//    0     Sol          -          0                     1.98E30     |
//    1     Mercurio	 0.2408	    47900                 3.28E23     |
//    2     Venus	     0.6152	    35000                 4.83E24     |
//    3     La Tierra	 1	        29800                 5.98E24     |
//    4     Marte	     1.8809	    24100                 6.40E23     |
//    5     Júpiter	     11.862	    13100                 1.90E27     |
//    6     Saturno	     29.458	    9600                  5.98E26     |
//    7     Urano	     84.01	    6800                  8.67E25     |
//    8     Neptuno	     164.79	    5400                  1.05E26     |
//    9     Plutón       248.54     4700                  1.30E22     |
//                                                                    |
//---------------------------------------------------------------------


#include <iostream>
#include <cmath>
#include <fstream>     // Para trabajar con ficheros

#define M_S 1.989E30   // Masa solar en kg
#define c 1.496E11     // Distancia entre la Tierra y el Sol
#define G 6.67384E-11  // Constante de gravitación universal
#define N 10           // Número de cuerpos (9 planetas + Sol)

using namespace std;

// FUNCIONES QUE SE VAN A UTILIZAR
double masa_solar(double m);
double distancia_UA(double d);
double velocidad_UA(double v);
int Newton(double masa[10], double (&r)[10][2], double (&a)[10][2]);
float Verlet(double masa[10], double (&r)[10][2], double (&v)[10][2], double (&a)[10][2], float h, float t);
float Tiempo_UA(float t);
void Periodo(string nombre_fichero, float periodo[N], int vueltas[N]);
double Energia(double masa[N], double r[N][2]);

/*---------------------------------------------------------------------
|                           FUNCIÓN PRINCIPAL                         |
---------------------------------------------------------------------*/
int main()
{
    // ------------------ DECLARACIÓN DE VARIABLES ---------------------
    
    float h;                    // Paso
    float t;                    // Tiempo
    float cambio_signo;         // Variable auxiliar para el cálculo de períodos
    float periodo[N];          // Periodo de cada planeta

    int i, j;                   // Contadores
    int vueltas[N];             // Número de vueltas que da cada planeta
    int elegir;                 // Para elegir las condiciones iniciales

    double r[10][2];            // Posición de los 9 planetas y el Sol
    double v[10][2];            // Velocidad de los 9 planetas y el Sol
    double a[10][2];            // Aceleración de los 9 planetas y el Sol
    double masa[10];            // Masa de los planetas y el Sol

    double matriz[N][5];       // Matriz donde guardo las condiciones iniciales
    double aux[N];              // Array auxiliar para el cálculo de periodos

    double energia;             // Energía potencial gravitatoria total

    ifstream inicial;           // Fichero con las condiciones iniciales
    ofstream fich_posicion;     // Fichero donde se van guardando las posiciones
    ofstream fich_energia;      // Fichero donde se guardan las energías totales

    // ----------------------- ABRIR FICHEROS --------------------------
    // Elegir si se quieren condiciones iniciales normales o con perturbaciones
    cout << "Elegir las condiciones iniciales (1=sin perturbar, 2=perturbadas)";
    cin >> elegir;
    if (elegir==1)
    {
        inicial.open("condiciones_iniciales.txt");
    }
    else if (elegir==2)
    {
        inicial.open("condiciones_iniciales_pert.txt");
    }
    else return 0;

    fich_posicion.open("posiciones_planetas.txt");
    fich_energia.open("energias_planetas.txt");

    // ------------------- CONDICIONES INICIALES -----------------------
    // Extraigo las condiciones iniciales de "condiciones_iniciales.txt"
    // Columna 0: masa
    // Columna 1: r_x (semieje mayor de la órbita)
    // Columna 2: r_y=0 
    // Columna 3: v_x=0
    // Columna 4: v_y (velocidad orbital)

    // Copio "condiciones_iniciales.txt" a un array "matriz[N][5]"
    if(!inicial.is_open()) cout << "Error al abrir el fichero";
    else
    {
        while (!inicial.eof())
        {
            for (i=0; i<N; i++)
            {  
                for (j=0; j<5; j++)
                {
                    inicial >> matriz[i][j]; // Copio el fichero a un array
                }
            }
        }
    }
    inicial.close();

    // Masas de los planetas y el Sol
    for(i=0; i<N; i++)
    {  
        masa[i]=matriz[i][0];           // Copio la masa (en kg) de la 1ª columna de "matriz"
        masa[i]=masa_solar(masa[i]);    // La escribo en masas solares
    }

    // Posición de los planetas y el Sol
    for(i=0; i<N; i++)
    {  
        r[i][0]=matriz[i][1];           // Componente x (en m)
        r[i][0]=distancia_UA(r[i][0]);  // La escribo en UA
        r[i][1]=matriz[i][2];           // Componente y es 0
        if (r[i][1] != 0) r[i][1]=distancia_UA(r[i][1]);  // La escribo en UA
    }
    
    // Velocidad de los planetas en km/s
    for(i=0; i<N; i++)
    {  
        v[i][0]=matriz[i][3];           // Componente x es 0
        if (v[i][0] != 0) v[i][0]=velocidad_UA(v[i][0]);  // La escribo en UA
        v[i][1]=matriz[i][4];           // Componente y (en m/s)
        v[i][1]=velocidad_UA(v[i][1]);  // La escribo en UA
    }
    
    // Definimos el paso e inicializamos el tiempo a t=0
    h = 0.01;
    t = 0.;

    // Inicializamos las variables utilizadas para el cálculo de periodos
    for(i=0; i<N; i++)
    {
        vueltas[i]=0;
        periodo[i]=0;
    }

    // -------------------CÁLCULO DE a(0)-----------------------
    Newton(masa, r, a);

    // --------------------- SIMULACIÓN ------------------------
    for(j=0; j<=156000; j++) // Simulo durante aprox. 248 años (periodo de Plutón)
    {
        // ---------- ENERGÍA POTENCIAL GRAVITATORIA -----------
        energia=Energia(masa, r);
        fich_energia << t << "\t" << energia << endl;

        // Guardar las posiciones para el vídeo y las energías
        if(j%100==0) // Quito algunas iteraciones para acortar el vídeo
        {
            for(i=0; i<N; i++)
            {
                fich_posicion << r[i][0] << ", " << r[i][1] << endl;
            }
            fich_posicion << endl;
        }
        // Guardo la componente "y" de la posición (para el período)
        for(i=1; i<N; i++)
        {
            aux[i]=r[i][1];
        }
        // Actualizo el tiempo y las posiciones/velocidades/aceleraciones
        t=Verlet(masa, r, v, a, h, t);

    // -------------------- CÁLCULO DE PERIODOS -------------------
        for(i=1;i<N;i++) // no cuento el Sol
        {
            cambio_signo=aux[i]*r[i][1];
            // Si r_y cambia de signo al actualizar t y r_x>0, sumo una vuelta
            if((cambio_signo<0)&&(r[i][0]>0))
            {
                vueltas[i]++;
                periodo[i]=t;
            }
        }

    } // Se acaba el bucle de iteraciones

    // Escribo los periodos en "periodos_planetas.txt"
    Periodo("periodos_planetas.txt", periodo, vueltas);
    
    // -------------------FIN DEL PROGRAMA----------------------
    // Cierro los ficheros
    fich_posicion.close();
    fich_energia.close();

    return 0;
}


/*---------------------------------------------------------------------
|                     FUNCIONES DE REESCALAMIENTO                     |
---------------------------------------------------------------------*/

// Función que transforma una masa en kg a unidades de masa solar
double masa_solar(double m)
{
    double masa;
    masa=m/M_S;
    return masa;
}

// Función que transforma una distancia en m a unidades astronómicas (UA)
double distancia_UA(double d)
{
    double distancia;
    distancia=d/c;
    return distancia;
}

// Función que transforma una velocidad en m/s a UA
double velocidad_UA(double v)
{
    double velocidad;
    velocidad=sqrt(c/(G*M_S))*v;
    return velocidad;
}


/*---------------------------------------------------------------------
|                          ALGORITMO DE VERLET                        |
---------------------------------------------------------------------*/

// Función que calcula la aceleración a partir de la 2ª LEY DE NEWTON
int Newton(double masa[10], double (&r)[10][2], double (&a)[10][2])
{
    int i, j;
    double resta_x, resta_y, modulo;

    // Inicializamos la aceleración a 0
    for(i=0; i<N; i++)
    {
        a[i][0]=0.;
        a[i][1]=0.;
    }
    // Calculamos la aceleración para cada planeta
    for(i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            if (i==j) continue; // Si i=j el denominador es 0, así que nos lo saltamos
            else
            {
                resta_x=r[i][0]-r[j][0];
                resta_y=r[i][1]-r[j][1];
                modulo=(resta_x*resta_x)+(resta_y*resta_y);
                a[i][0]=a[i][0]-(masa[j]*resta_x/pow(modulo,1.5));     // Componente x de la aceleración
                a[i][1]=a[i][1]-(masa[j]*resta_y/pow(modulo,1.5));     // Componente y de la aceleración
            }
        }
    }
    return 0;
}

// Función que itera el ALGORITMO DE VERLET a partir de r(0), v(0) y a(0)
float Verlet(double masa[10], double (&r)[10][2], double (&v)[10][2], double (&a)[10][2], float h, float t)
{ 
    double w[10][2];    // Vector auxiliar "w" del algoritmo de Verlet
    int i;              // Contador

    // Paso 1: calcular w(t+h)=v(t)+(0.5*h)a(t)
    for(i=0; i<N; i++)
    {
        w[i][0]=v[i][0]+0.5*h*a[i][0];
        w[i][1]=v[i][1]+0.5*h*a[i][1];
    }
    // Paso 2: calcular r(t+h)=r(h)+hv(h)+(0.5*h*h)a(h)
    for(i=0; i<N; i++)
    {
        r[i][0]=r[i][0]+h*w[i][0];
        r[i][1]=r[i][1]+h*w[i][1];
    }
    // Paso 3: calcular a(t+h) con Newton usando los nuevos valores de r(t+h)
    Newton(masa, r, a);
    // Paso 4: calcular v(t+h)=w(t+h)+(0.5*h)a(t+h) con los nuevos valores de a(t+h)
    for(i=0; i<N; i++)
    {
        v[i][0]=w[i][0]+0.5*h*a[i][0];
        v[i][1]=w[i][1]+0.5*h*a[i][1];
    }
    // Paso 5: devolver el valor de t+h
    return t+h; 
}

/*---------------------------------------------------------------------
|                FUNCIONES PARA EL CÁLCULO DE PERIODOS                |
---------------------------------------------------------------------*/

// Función que transforma un tiempo en UA a años
float Tiempo_UA(float t)
{
    float tiempo;
    tiempo=(sqrt(pow(c,3)/(G*M_S))/(3600*24*365))*t;
    return tiempo;
}

// Función que escribe el PERIODO de los planetas en un fichero "periodos_planetas.txt"
/* Tengo que dividir el tiempo guardado en periodo[N] entre el número de vueltas
guardado en vueltas[N]. Luego paso ese tiempo a segundos usando Tiempo_UA(t) */
void Periodo(string nombre_fichero, float periodo[N], int vueltas[N])
{
    int i;
    ofstream fich_periodos;       // Fichero donde se van guardando los periodos
    double periodo_planeta[N];    // Valor de los periodos en segundos
    float error_relativo[N];      // Error relativo de los periodos

    fich_periodos.open(nombre_fichero);

    for (i=1; i<N; i++)
    {
        // Calculo el periodo
        periodo_planeta[i]=periodo[i]/vueltas[i];
    }

        // Calculo el error relativo en porcentaje (%)
        error_relativo[1]=100*(Tiempo_UA(periodo_planeta[1])-0.24)/0.24;      // Mercurio
        error_relativo[2]=100*(Tiempo_UA(periodo_planeta[2])-0.615)/0.615;    // Venus
        error_relativo[3]=100*(Tiempo_UA(periodo_planeta[3])-1)/1;            // Tierra
        error_relativo[4]=100*(Tiempo_UA(periodo_planeta[4])-1.88)/1.88;      // Marte
        error_relativo[5]=100*(Tiempo_UA(periodo_planeta[5])-11.86)/11.86;    // Júpiter
        error_relativo[6]=100*(Tiempo_UA(periodo_planeta[6])-29.46)/29.46;    // Saturno
        error_relativo[7]=100*(Tiempo_UA(periodo_planeta[7])-84.01)/84.01;    // Urano
        error_relativo[8]=100*(Tiempo_UA(periodo_planeta[8])-164.79)/164.79;  // Neptuno
        error_relativo[9]=100*(Tiempo_UA(periodo_planeta[9])-247.92)/247.92;  // Plutón

    for (i=1; i<N; i++)
    {
        // Guardo en un fichero el valor del periodo en años
        fich_periodos << "Planeta " << i << ": " << Tiempo_UA(periodo_planeta[i]) << " años / Error relativo: "<< std::abs(error_relativo[i]) << " %" << endl;
    }

    fich_periodos.close();
    return;
}


/*---------------------------------------------------------------------
|                FUNCIONES PARA EL CÁLCULO DE ENERGÍAS                |
---------------------------------------------------------------------*/

// Función que devuelve la ENERGÍA potencial gravitatoria total y la
// escribe en un fichero "energias_planetas.txt"
double Energia(double masa[N], double r[N][2])
{
    double modulo, resta_x, resta_y;
    double energia[N];     // Guarda la energía gravitatoria de cada planeta
    double energia_total;  // Guarda la suma de todas las energías
    int i, j;

    // Inicializo energia[N] y energia_total
    for(i=0; i<N; i++)
    {
        energia[i]=0;
    }
    energia_total=0;

    for(i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            if (i==j) continue; // Si i=j el denominador es 0, así que nos lo saltamos
            else
            {
                resta_x=r[i][0]-r[j][0];
                resta_y=r[i][1]-r[j][1];
                modulo=sqrt((resta_x*resta_x)+(resta_y*resta_y));
                energia[i]=energia[i]+(masa[j]/modulo);
            }
        }
        energia_total=energia_total-(masa[i]*energia[i]);
    }

    return energia_total;
}