#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biblioteca.h" //Biblioteca com funções necessárias

/*
    Para executar o código em sistemas windows você deve:
    - Abrir a pasta com o arquivo executável (.exe). Na pasta original:IC.2020.10.1 -> bin -> Debug;
    - Clicar em "Arquivo" no canto superior esquerdo e selecionar "Abrir o Windows PowerShell";
    - Com o PowerShell aberto, digite .\(nome do arquivo). No arquivo original: .\IC.2020.10.exe
    - Com a chamada do programa, digitar novamente .\(nome do arquivo) seguido das variáveis de lançamento pedidas separadas por espaço.

*/

//Código para cálculo da trajetória e velocidade em coordenadas esféricas do voo de um balão meteorológico com hélio
int main(int argc, char* argv[])
{
    double r, theta, lambda, vol, m, alt, h, r_b, r_b0, e, e_0, r_bf, m_h, m_s, g, T, T0, P, P_h, ro, ro_h0, ro_h, r_ant;
    double d_r, d_lambda, d_theta, vvento_r, vvento_theta, vvento_lambda, drb_dpm, beta_e;
    double t = 0.0; //Tempo
    int cont_p = 0.0;//Contador para impressão do resultado de 1000 a 1000 passos
    h = 0.01;//Passo
    int result;

    if(argc != 8){
        printf("Usar: IC.2020.10.27.exe altitude(m) latitude(graus) longitude(graus) massa_balao(balão + payload [kg]) empuxo(kgf) Altitude de flutuação(m) Diâmetro da válvula (m)\n\n");
        exit(1);
    }

    //Aplicação das variáveis dadas pelo usuário no código
    double alt_0 = atof(argv[1]); //altitude de lançamento. Ex: 850.0;
    double latitude = atof(argv[2]);//latitude de lançamento. Ex: -22.00548;
    double longitude = atof(argv[3]);//longitude de lançamento. Ex: -47.93405;
    m = atof(argv[4]); //massa balão + payload. Ex: 6.2615;
    double empuxo = atof(argv[5]); //Ex: 8.390508;
    double alt_flut = atof(argv[6]); //Ex: 25000
    double d = atof(argv[7]); //Diâmetro da válvula. Ex. 0.02

    //Transformação dos dados para variávies do programa
    r = alt_0 +R_T;
    theta = M_PI/2 - latitude*M_PI/180;
    lambda = longitude*M_PI/180;

    //Determina-se derivadas iniciais
    d_r = 0.0;
    d_theta = 0.0;
    d_lambda = 0.0;

    //Determina-se as velocidades iniciais do vento (cte)
    vvento_r = 0.0;
    vvento_theta = 0.0;
    vvento_lambda = 0.0;


    FILE *arq; // chamando arquivo para escrita de dados
    arq = fopen("Teste1.txt","w"); //Abrindo arquivo para escrita de dados

    //Impressão do nome das variáveis no início do arquivo
    result = fprintf(arq,"#       Tempo      Altitude       Pressão   Empuxo    Volume     Densidade\n");

    //Chamada de funções e determinação de variávies
    g = grav(r); //aceleração da gravidade local
    T0 = T= 22.0+273.15; //Alteracao 2020-11-05 Temp(r); //temperatura local
    P =  94360; //Alteracao 2020-11-05 Pres(r, T, g); //pressão local
    ro = P/(Ra*T); //densidade atmosférica local
    ro_h = P/(R_h*T); //densidade local do hélio
    vol = empuxo/ro; //volume inicial do balão. Ex: 7.5464;
    m_h = ro_h * vol; //massa de hélio
    m_s = m - m_h;
    e_0 = e = 0.0005; //espessura da membrana de látex. Deve ser alterado com valores medidos!!
    r_b = r_b0 = raio_bi(ro_h, m_h);


    // Determinação das "inclinações" dos intervalos para o método de Hunge Kutta através de vetores de 6 elementos
    double k1[8],k2[8],k3[8],k4[8];

    empuxo = ro*vol*(g+(r*OMEGA*OMEGA*sin(theta)*sin(theta)));
    double peso = m*g;
    //Impressão dos dados iniciais no arquivo
    result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T, P, empuxo,  vol, ro);

    //Cálculo dos vetores k's para cada elemento/variável posição, velocidade e raio do balão
    while (r - R_T <= 32000){

        //cálculo das variáves atmosféricas para cada altitude
        alt = r - R_T;
        g = grav(r);
        update_atm(alt_0,T0,r,r_ant,g,&T,&P,&ro);
        r_ant = r;
        m_h = m - m_s;
        ro_h = m_h/vol; //densidade do hélio
        P_h = ro_h*R_h*T; //pressão dentro do balão adotando temperatura interna = externa

        //Volume do balão - aproximação esférica
        vol = 4.0*M_PI*pow(r_b,3.0)/3.0;


        //Modelo teórico de elasticidade da membrana de látex
        e = e_0*r_b0*r_b0/(r_b*r_b); //cálculo da redução da espessura da membrana conforme crescimento do volume/raio do balão
        beta_e = 50990.23198*pow(r_b, -4.885887); //coeficiente de elasticidade X espessura
        drb_dpm = r_b*r_b/(2.0*beta_e - 3.0*r_b*(P_h-P)); //derivada do raio do balão pela pressão manométrica

        //calculo de k1 para cada uma das variáveis
        k1[0] = h*drb_dt(r, d_r, T, P, r_b, ro_h, g, m, P_h, drb_dpm, alt_0, T0);
        k1[1] = h*d_r;
        k1[2] = h*d_theta;
        k1[3] = h*d_lambda;
        k1[4] = h*d2r(r,theta,d_r,d_theta, d_lambda, m, vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[5] = h*d2theta(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[6] = h*d2lambda(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[7] = h*fluxo_massa(ro_h, d, P, P_h);

        //calculo k2
        k2[0] = h*drb_dt(r+k1[1]/2.0, d_r+k1[4]/2.0, T, P, r_b+k1[0]/2.0, ro_h, g, m+k1[7]/2.0, P_h, drb_dpm, alt_0, T0);
        k2[1] = h*(d_r+k1[4]/2.0);
        k2[2] = h*(d_theta+k1[5]/2.0);
        k2[3] = h*(d_lambda+k1[6]/2.0);
        k2[4] = h*d2r(r+k1[1]/2.0,theta+k1[2]/2.0,d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m+k1[7]/2.0,vol, g, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[5] = h*d2theta(r+k1[1]/2.0,theta+k1[2]/2.0,d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m+k1[7]/2.0,vol, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[6] = h*d2lambda(r+k1[1]/2.0, theta+k1[2]/2.0, d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m+k1[7]/2.0,vol, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[7] = h*fluxo_massa(ro_h, d, P, P_h);

        //calculo de k3
        k3[0] = h*drb_dt(r+k2[1]/2.0, d_r+k2[4]/2.0, T, P, r_b+k2[0]/2.0, ro_h, g, m+k2[7]/2.0, P_h, drb_dpm, alt_0, T0);
        k3[1] = h*(d_r+k2[4]/2.0);
        k3[2] = h*(d_theta+k2[5]/2.0);
        k3[3] = h*(d_lambda+k2[6]/2.0);
        k3[4] = h*d2r(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m+k2[7]/2.0,vol, g, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[5] = h*d2theta(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m+k2[7]/2.0,vol, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[6] = h*d2lambda(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m+k2[7]/2.0,vol, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[7] = h*fluxo_massa(ro_h, d, P, P_h);

        //calculo de k4
        k4[0] = h*drb_dt(r+k3[1], d_r+k3[4], T, P, r_b+k3[0], ro_h, g, m+k3[7], P_h, drb_dpm, alt_0, T0);
        k4[1] = h*(d_r+k3[4]);
        k4[2] = h*(d_theta+k3[5]);
        k4[3] = h*(d_lambda+k3[6]);
        k4[4] = h*d2r(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m+k3[7],vol, g, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[5] = h*d2theta(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m+k3[7],vol, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[6] = h*d2lambda(r+k3[1],theta+k3[2], d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m+k3[7],vol, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);

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

        //Vazão
        if ((alt>alt_flut-1000)&&(d_r>0)){
            m -= (k1[7]+2.0*k2[7]+2.0*k3[7]+k4[7])/6.0;
            alt = r - R_T;
            printf("%13.8lf %13.8lf\n", alt, m);
        continue;
        }

        empuxo = ro*vol*(g+(r*OMEGA*OMEGA*sin(theta)*sin(theta)));
        peso = m*g;

        //Impressão das variáveis no arquivo
        if (cont_p==2){ // A cada 1000 loops do while os valores das variáveis são impressas no arquivo
            if(arq==NULL){ //Mensagem caso haja erro na abertura do arquivo
                printf("Problemas na abertura do arquivo\n");
                system("pause");
                exit(1);
            }
            printf("%13.8lf %13.8lf %13.8lf %13.8lf\n",t , r, d_r, m);//impressão de dados desejados na tela
            //impressão de dados no arquivo
            result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T, P, empuxo,  vol,ro);
            if(result<0){
                printf("Erro na escrita do arquivo\n");
            }
            cont_p = 0; //Zera o contador
        }
        }
        fclose(arq); //Fecha o arquivo
        system("pause");
    return 0;
}

