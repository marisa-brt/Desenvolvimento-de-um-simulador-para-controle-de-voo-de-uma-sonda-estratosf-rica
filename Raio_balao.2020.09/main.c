#define Ca 0.34037379532479679571
#define Rh 2077.1 //Constante individual do hélioderivada do raio do balão pela pressão interna
#define Ra 287.54 // Constante universal dos gases tomando a massa molar do ar constante
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Determinação do raio do balão conforme densidade atmosférica (tomando v_r cte)

int main(){
double d_r, g, m, r_b, ro;
double a, b, r_bf, n, vol_b;
double alt, T, T0, P, P0, grad_t, C, P_h, m_h, ro_h, dt, alt_f, result;

double lamb, drdp, esp_0, vol_m0, beta_e;

//beta_e = beta * espessura

FILE *arq; // chamando arquivo para escrita de dados
arq = fopen("r_b_V3.dat","w"); //Abrindo arquivo para escrita de dados
result = fprintf(arq,"    Altitude       Raio_b        P_h - P       dR/dp       Beta_e\n");

esp_0 = 0.0005; //espessura inicial do balão

//Dados G3
d_r = 4.873; //velocidade do balão - mantida constante na subida
g = 9.8168; //aceleração da gravidade
m = 6.2615; //massa do balão

T0 = 295.15; //temperatura no lançamento
P0 = 94360; //pressão atmosférica no lançamento

ro = P0/Ra/T0; //densidade do ar
ro_h = P0/Rh/T0; //densidade do hélio (nas mesmas condições de temperatura e pressão)

vol_b = 7.54644637536312716197; //volume inicial do balão
r_b = pow(3.0*vol_b/4.0/M_PI,1.0/3.0); //raio inicial do balão

m_h = ro_h*vol_b; //massa de hélio

alt = 850.0; //altitude de lançamento
dt = 1.0;
a = Ca*3.0*d_r*d_r/(8.0*g);
b = 3.0*m/(4.0*M_PI);

while(alt <= 11000){
        n = 0;
        //método de Newton para determinar o raio do balão dada a densidade do ar
        while (n <= 10){
            r_bf = r_b - (pow(r_b,3.0)- a*r_b*r_b -(b/ro))/(3.0*r_b*r_b - 2.0*a*r_b);
            r_b = r_bf;
            n +=1;
        }

        //T = -0.0049687*alt + 302.91868;
        //grad_t = -6.5e-3; //gradiente térmico específico da camada atmosférica em K/m (ISA)

        grad_t = -4.9687e-3; //medido no voo
        T = T0 + grad_t*(alt-850.0); //utilizando dT/dz medido no balao.
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
        ro = P/(Ra*T);
        ro_h = m_h/(4.0*M_PI*pow(r_b,3.0)/3.0);
        P_h = ro_h*Rh*T;
        alt_f = alt + d_r*dt;
        alt = alt_f;

        //Modelo termodinâmico para determinação da derivada do raio do balão pela pressão interna
        lamb = (r_b - a)/(3.0*r_b - 2.0*a)*(C-1.0);
        drdp = (r_b*lamb)/(C*P - P_h*(1.0 + 3.0*lamb));//derivada do raio do balão pela pressão interna
        //e = e_0*r_b*r_b/(r_b0*r_b0);//determinação da espessura
        beta_e = r_b/2.0*(r_b/drdp + 3.0*(P_h-P));

        printf("%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",alt, r_b, P_h-P, drdp, beta_e);
        result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",alt, r_b, P_h-P, drdp, beta_e);
}


    }



