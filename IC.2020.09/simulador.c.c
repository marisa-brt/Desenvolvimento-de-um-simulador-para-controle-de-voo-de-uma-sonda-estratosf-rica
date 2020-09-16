#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biblioteca.h" //Biblioteca com funções necessárias

//Código para cálculo da trajetória e velocidade em coordenadas esféricas do voo de um balão meteorológico com hélio
int main(int argc, char* argv[])
{
    double r, theta, lambda, vol, m, alt, h, r_b, r_b0, e, e_0, r_bf, m_h, g, T, P, P_h, ro, ro_h0, ro_h;
    double d_r, d_lambda, d_theta, vvento_r, vvento_theta, vvento_lambda, drb_dpm, beta_e;
    double t = 0.0; //Tempo
    int cont_p = 0.0;//Contador para impressão do resultado de 1000 a 1000 passos
    h = 0.01;//Passo
    int result;

    //if(argc != 6){
        //printf("Usar: IC.2020.08.exe altitude latitude(graus) longitude(graus) massa_balao(balão + payload [kg]) empuxo(kgf)\n\n");
        //exit(1);
    //}

    //Aplicação das variáveis dadas pelo usuário no código
    alt = 850.0;//atof(argv[1]);
    double latitude = -22.00548;//atof(argv[2]);
    double longitude = -47.93405;//atof(argv[3]);
    m = 6.2615;//atof(argv[4]);
    //double empuxo = atof(argv[5]) + m;

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
    arq = fopen("Teste4.txt","w"); //Abrindo arquivo para escrita de dados
    //printf("Defina as coordenadas (theta e lambda), o volume do balão,a massa e a altitude de lançamento \n");
    //scanf("%lf %lf %lf %lf %lf", &theta, &lambda, &vol, &m, &alt);//Informações definidas pelo usuário

    result = fprintf(arq,"#       Tempo      Altitude       Vel rad     Densid    Volume    Raio balão\n");

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
    vol = 7.54644637536312716197;//empuxo/ro; //volume inicial do balão
    m_h = ro_h * vol; //massa de hélio
    e_0 = e = 0.0005;
    r_b = r_b0 = raio_bi(ro_h, m_h);

    /*
    //Leitura de valores do raio do balão e do beta em arquivo
    FILE *arq2 = fopen("rbteste.dat","rb");
    double rbarq, betaarq, r_ba[2084], betaa[2084];
    int cont_p2 = 1;

        if(arq2 == NULL){
            printf("Arquivo de raio de balão x beta não abre!");
        }
        else{
            while((fscanf(arq2, "   %lf %lf", &rbarq, &betaarq)) != EOF){
                r_ba[cont_p2] = rbarq;
                betaa[cont_p2] = betaarq;
                cont_p2 += 1;
            }
        }
        int i = 1;
        for(; i<cont_p2; ++i){
            //printf(" %13.8lf %13.8lf\n", r_ba[i], betaa[i]);
        }
        fclose(arq2);
        i = 1;*/

    // Determinação das "inclinações" dos intervalos para o método de Hunge Kutta através de vetores de 6 elementos
    double k1[7],k2[7],k3[7],k4[7];
    //Cálculo dos vetores k's para cada elemento/variável posição e velocidade

    result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T,  d_r,  ro, vol,  r_b);//troquei d_theta e d_lambda por ro e vol

    while ( r - R_T <= 32000){
        g = grav(r);
        T = Temp(r);
        P = Pres(r, T, g);
        ro = P/(Ra*T);
        vol = 4.0*M_PI*pow(r_b,3.0)/3.0;//m_h/ro_h;
        ro_h = m_h/vol;
        P_h = ro_h*R_h*T;
        e = e_0*r_b0*r_b0/(r_b*r_b);

        /*
        //Determinação de beta
        while (fabs(r_b - r_ba[i])<=0.0001){
            if (fabs(r_b - r_ba[i-1]<fabs(r_b - r_ba[i]))){
                beta = betaa[i-1];
                break;
                }
            else beta = betaa[i];
            i++;
        }*/
        beta_e = 50990.23198*pow(r_b, -4.885887);
        drb_dpm = r_b*r_b/(2.0*beta_e - 3.0*r_b*(P_h-P));

        //calculo de k1 para cada uma das variáveis
        k1[0] = h*drb_dt(r, d_r, T, P, r_b, ro_h, g, m, P_h, drb_dpm);
        k1[1] = h*d_r;
        k1[2] = h*d_theta;
        k1[3] = h*d_lambda;
        k1[4] = h*d2r(r,theta,d_r,d_theta, d_lambda, m, vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[5] = h*d2theta(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[6] = h*d2lambda(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);

        //calculo k2
        k2[0] = h*drb_dt(r+k1[1]/2.0, d_r+k1[4]/2.0, T, P, r_b+k1[0]/2.0, ro_h, g, m, P_h, drb_dpm);
        k2[1] = h*(d_r+k1[4]/2.0);
        k2[2] = h*(d_theta+k1[5]/2.0);
        k2[3] = h*(d_lambda+k1[6]/2.0);
        k2[4] = h*d2r(r+k1[1]/2.0,theta+k1[2]/2.0,d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m,vol, g, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[5] = h*d2theta(r+k1[1]/2.0,theta+k1[2]/2.0,d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m,vol, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[6] = h*d2lambda(r+k1[1]/2.0, theta+k1[2]/2.0, d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m,vol, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);

        //calculo de k3
        k3[0] = h*drb_dt(r+k2[1]/2.0, d_r+k2[4]/2.0, T, P, r_b+k2[0]/2.0, ro_h, g, m, P_h, drb_dpm);
        k3[1] = h*(d_r+k2[4]/2.0);
        k3[2] = h*(d_theta+k2[5]/2.0);
        k3[3] = h*(d_lambda+k2[6]/2.0);
        k3[4] = h*d2r(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m,vol, g, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[5] = h*d2theta(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m,vol, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[6] = h*d2lambda(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m,vol, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);

        //calculo de k4
        k4[0] = h*drb_dt(r+k3[1], d_r+k3[4], T, P, r_b+k3[0], ro_h, g, m, P_h, drb_dpm);
        k4[1] = h*(d_r+k3[4]);
        k4[2] = h*(d_theta+k3[5]);
        k4[3] = h*(d_lambda+k3[6]);
        k4[4] = h*d2r(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m,vol, g, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[5] = h*d2theta(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m,vol, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[6] = h*d2lambda(r+k3[1],theta+k3[2], d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m,vol, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);

        t += h;// O tempo é acrescido pelo passo h a cada loop
        cont_p++; //Contador de impressão se soma

        //Cálculo das variáveis
        r_b += (k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0])/6.0;
        r += (k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])/6.0;
        theta += (k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2])/6.0;
        lambda += (k1[3]+2.0*k2[3]+2.0*k3[3]+k4[3])/6.0;
        d_r += (k1[4]+2.0*k2[4]+2.0*k3[4]+k4[4])/6.0;
        d_theta += (k1[5]+2.0*k2[5]+2.0*k3[5]+k4[5])/6.0;
        d_lambda += (k1[6]+2.0*k2[6]+2.0*k3[6]+k4[6])/6.0;

        printf("%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",t , r, d_r, d_theta, d_lambda);

        if (cont_p==10){ // A cada 100 loops do while os valores das variáveis são impressas no arquivo
            if(arq==NULL){ //Mensagem caso haja erro na abertura do arquivo
                printf("Problemas na abertura do arquivo\n");
                system("pause");
                exit(1);
            }
            //impressão de dados no arquivo
            result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T,  d_r,  ro, vol,  r_b);//troquei d_theta*r,  d_lambda*r*sin(theta) por ro e vol
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
