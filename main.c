#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biblioteca.h" //Biblioteca com funções necessárias

//Código para cálculo da trajetória e velocidade em coordenadas esféricas do voo de um balão meteorológico com hélio
int main(int argc, char* argv[])
{
    double r, theta,lambda,vol, m, alt, h, r_b, m_h, g, T, P, ro, ro_h0, ro_h;
    double d_r, d_lambda, d_theta, vvento_r, vvento_theta, vvento_lambda;
    double t = 0.0; //Tempo
    int cont_p = 0.0;//Contador para impressão do resultado de 1000 a 1000 passos
    h = 0.01;//Passo
    int result;

    if(argc != 6){
        printf("Usar: IC.2020.exe altitude latitude(graus) longitude(graus) massa_balao(balão + payload [kg]) empuxo(kgf)\n\n");
        exit(1);
    }

    //Aplicação das variáveis dadas pelo usuário no código
    alt = atof(argv[1]);
    double latitude = atof(argv[2]);
    double longitude = atof(argv[3]);
    m = atof(argv[4]);
    double empuxo = atof(argv[5]);

    r = alt+R_T;
    theta = M_PI/2 - latitude*M_PI/180;
    lambda = longitude*M_PI/180;

    //Determina-se derivadas iniciais nulas
    d_r = 0.0;
    d_theta = 0.0;
    d_lambda = 0.0;

    //Determina-se as velocidade iniciais do vento (cte)
    vvento_r = 0.0;
    vvento_theta = 0.0;
    vvento_lambda = 0.0;


    FILE *arq; // chamando arquivo para escrita de dados
    arq = fopen("Teste.txt","w"); //Abrindo arquivo para escrita de dados
    //printf("Defina as coordenadas (theta e lambda), o volume do balão,a massa e a altitude de lançamento \n");
    //scanf("%lf %lf %lf %lf %lf", &theta, &lambda, &vol, &m, &alt);//Informações definidas pelo usuário

    result = fprintf(arq,"#       Tempo      Altitude       Vel rad     Vel theta    Vel lambda    Raio balão\n");

    //theta = 1.54;//-0.3839724;
    //lambda = 1.54;//2.3387412;
    //vol_0 = 5/1.1226;
    //m = 6.26;
    //alt = 852.0;
    //vol = 7.666;
    //ro_h0 = 0.151;
    //g = 9.8168;
    //ro = 1.094;
    g = grav(r); //aceleração da gravidade local
    T = Temp(r); //temperatura local
    P = Pres(r, T, g); //pressão local
    ro = P/(Ra*T); //densidade atmosférica local
    ro_h = P/(R_h*T); //densidade local do hélio
    vol = empuxo/ro; //volume inicial do balão
    m_h = ro_h * vol; //massa de hélio
    r_b = raio_bi(ro_h, m_h);//r_b=r_bi


    // Determinação das "inclinações" dos intervalos para o método de Hunge Kutta através de vetores de 6 elementos
    double k1[7],k2[7],k3[7],k4[7];
    //Cálculo dos vetores k's para cada elemento/variável posição e velocidade

    result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T,  d_r,  d_theta,  d_lambda,  r_b);

    while ( r - R_T <= 32000){
        //g = grav(r);
        T = Temp(r);
        P = Pres(r, T, g);
        //ro = P/(Ra*T);
        //ro_h = P/(R_h*T);
        //vol = m_h/ro_h;
        //r_b = pow((3.0*vol/(4.0*M_PI)),(1.0/3.0));
        //r_b = raio_b(P, T, m_h, ro_h0, r_b);

        //calculo de k1 para cada uma das variáveis
        k1[1] = h*d_r;
        k1[2] = h*d_theta;
        k1[3] = h*d_lambda;
        k1[4] = h*d2r(r,theta,d_r,d_theta, d_lambda, m, vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[5] = h*d2theta(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[6] = h*d2lambda(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);

        //calculo k2
        k2[1] = h*(d_r+k1[4]/2.0);
        k2[2] = h*(d_theta+k1[5]/2);
        k2[3] = h*(d_lambda+k1[6]/2);
        k2[4] = h*d2r(r+k1[1]/2,theta+k1[2]/2,d_r+k1[4]/2,d_theta+k1[5]/2, d_lambda+k1[6]/2,m,vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k2[5] = h*d2theta(r+k1[1]/2,theta+k1[2]/2,d_r+k1[4]/2,d_theta+k1[5]/2, d_lambda+k1[6]/2,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k2[6] = h*d2lambda(r+k1[1]/2, theta+k1[2]/2, d_r+k1[4]/2,d_theta+k1[5]/2, d_lambda+k1[6]/2,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);

        //calculo de k3
        k3[1] = h*(d_r+k2[4]/2.0);
        k3[2] = h*(d_theta+k2[5]/2);
        k3[3] = h*(d_lambda+k2[6]/2);
        k3[4] = h*d2r(r+k2[1]/2,theta+k2[2]/2,d_r+k2[4]/2,d_theta+k2[5]/2, d_lambda+k2[6]/2,m,vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k3[5] = h*d2theta(r+k2[1]/2,theta+k2[2]/2,d_r+k2[4]/2,d_theta+k2[5]/2, d_lambda+k2[6]/2,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k3[6] = h*d2lambda(r+k2[1]/2,theta+k2[2]/2,d_r+k2[4]/2,d_theta+k2[5]/2, d_lambda+k2[6]/2,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);

        //calculo de k4
        k4[1] = h*(d_r+k3[4]);
        k4[2] = h*(d_theta+k3[5]);
        k4[3] = h*(d_lambda+k3[6]);
        k4[4] = h*d2r(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m,vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k4[5] = h*d2theta(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k4[6] = h*d2lambda(r+k3[1],theta+k3[2], d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);

        t += h;// O tempo é acrescido pelo passo h a cada loop
        cont_p++; //Contador de impressão se soma

        //Cálculo das variáveis
        r += (k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        theta += (k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        lambda += (k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
        d_r += (k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;
        d_theta += (k1[5]+2*k2[5]+2*k3[5]+k4[5])/6;
        d_lambda += (k1[6]+2*k2[6]+2*k3[6]+k4[6])/6;

        printf("%13.8lf  %13.8lf %13.8lf %13.8lf %13.8lf\n",t , r, d_r, d_theta, d_lambda);

        if (cont_p==100){ // A cada 100 loops do while os valores das variáveis são impressas no arquivo
            if(arq==NULL){ //Mensagem caso haja erro na abertura do arquivo
                printf("Problemas na abertura do arquivo\n");
                system("pause");
                exit(1);
            }
            //impressão de dados no arquivo
            result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T,  d_r,  d_theta*r,  d_lambda*r*sin(theta),  r_b);
            if(result<0){
                printf("Erro na escrita do arquivo\n");
            }//
            cont_p = 0; //Zera o contador
        }
        }
        fclose(arq); //Fecha o arquivo
        system("pause");
    return 0;
}
