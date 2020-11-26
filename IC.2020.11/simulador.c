#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biblioteca.h" //Biblioteca com fun��es necess�rias

/*
    Para executar o c�digo em sistemas windows voc� deve:
    - Abrir a pasta com o arquivo execut�vel (.exe). Na pasta original:IC.2020.10.1 -> bin -> Debug;
    - Clicar em "Arquivo" no canto superior esquerdo e selecionar "Abrir o Windows PowerShell";
    - Com o PowerShell aberto, digite .\(nome do arquivo). No arquivo original: .\IC.2020.10.exe
    - Com a chamada do programa, digitar novamente .\(nome do arquivo) seguido das vari�veis de lan�amento pedidas separadas por espa�o.

*/

//C�digo para c�lculo da trajet�ria e velocidade em coordenadas esf�ricas do voo de um bal�o meteorol�gico com h�lio
int main(int argc, char* argv[])
{
    double r, theta, lambda, vol, m, alt, h, r_b, r_b0, e, e_0, m_h, g, T, T0, P, P_h, ro, ro_h, r_ant;
    double d_r, d_lambda, d_theta, vvento_r, vvento_theta, vvento_lambda, drb_dpm, beta_e;
    double Q, Cap, T_h;
    double t = 0.0; //Tempo
    int cont_p = 0.0;//Contador para impress�o do resultado de 1000 a 1000 passos
    h = 0.01;//Passo
    int result;

    /*if(argc != 9){
        printf("Usar: IC.2020.11.exe altitude latitude(graus) longitude(graus) massa_balao[kg] carga_paga(payload[kg]) empuxo(kgf) Altitude de flutua��o(m) Di�metro da v�lvula (m)\n\n");
        exit(1);
    }*/

    //Aplica��o das vari�veis dadas pelo usu�rio no c�digo
    double alt_0 = atof(argv[1]); //altitude de lan�amento. Ex: 850.0;
    double latitude = atof(argv[2]);//latitude de lan�amento. Ex: -22.00548;
    double longitude = atof(argv[3]);//longitude de lan�amento. Ex: -47.93405;
    double m_b = atof(argv[4]); //massa bal�o(sem h�lio). Ex: 1.6;
    double m_s = atof(argv[5]); //payload. Ex: 3.5115;
    double empuxo = atof(argv[6]); //Ex: 8.390508;
    double alt_flut = atof(argv[6]); //Ex: 25000
    double d = atof(argv[7]); //Di�metro da v�lvula. Ex. 0.02

    //Transforma��o dos dados para vari�vies do programa
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
    arq = fopen("Teste2.txt","w"); //Abrindo arquivo para escrita de dados

    //Impress�o do nome das vari�veis no in�cio do arquivo
    result = fprintf(arq,"#       Tempo      Altitude    Velocidade    Massa         Temperatura        Massa h�lio     T_h\n");

    //Chamada de fun��es e determina��o de vari�vies
    g = grav(r); //acelera��o da gravidade local
    T0 = T = 22.0+273.15; //Alteracao 2020-11-05 Temp(r); //temperatura local
    T_h = T + 0.05; //Temperatura do h�lio
    P =  94360; //Alteracao 2020-11-05 Pres(r, T, g); //press�o local
    ro = P/(Ra*T); //densidade atmosf�rica local
    ro_h = P/(R_h*T); //densidade local do h�lio
    vol = empuxo/ro; //volume inicial do bal�o. Ex: 7.5464;
    m_h = ro_h * vol; //massa de h�lio
    m = m_b + m_s + m_h;
    e_0 = e = 0.0005; //espessura da membrana de l�tex. Deve ser alterado com valores medidos!!
    r_b = r_b0 = raio_bi(ro_h, m_h);
    Q = 0.0;


    // Determina��o das "inclina��es" dos intervalos para o m�todo de Hunge Kutta atrav�s de vetores de 6 elementos
    double k1[9],k2[9],k3[9],k4[9];

    empuxo = ro*vol*(g+(r*OMEGA*OMEGA*sin(theta)*sin(theta)));
    double peso = m*g;
    //Impress�o dos dados iniciais no arquivo
    result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T, d_r, m,  T, m_h, T_h);

    //C�lculo dos vetores k's para cada elemento/vari�vel posi��o, velocidade e raio do bal�o
    while (r - R_T <= 32000){
        //c�lculo das vari�ves atmosf�ricas para cada altitude
        alt = r - R_T;
        g = grav(r);
        update_atm(alt_0,T0,r,r_ant,g,&T,&P,&ro);
        r_ant = r;
        m_h = m - m_s - m_b;
        ro_h = m_h/vol; //densidade do h�lio
        P_h = ro_h*R_h*T_h; //press�o dentro do bal�o adotando temperatura interna = externa

        //Volume do bal�o - aproxima��o esf�rica
        vol = 4.0*M_PI*pow(r_b,3.0)/3.0;

        //Modelo te�rico de elasticidade da membrana de l�tex
        e = e_0*r_b0*r_b0/(r_b*r_b); //c�lculo da redu��o da espessura da membrana conforme crescimento do volume/raio do bal�o
        beta_e = 50990.23198*pow(r_b, -4.885887); //coeficiente de elasticidade X espessura
        drb_dpm = r_b*r_b/(2.0*beta_e - 3.0*r_b*(P_h-P)); //derivada do raio do bal�o pela press�o manom�trica

        //calculo de k1 para cada uma das vari�veis
        k1[0] = h*drb_dt(r, d_r, T, P, r_b, ro_h, g, m, P_h, drb_dpm, alt_0, T0);
        k1[1] = h*d_r;
        k1[2] = h*d_theta;
        k1[3] = h*d_lambda;
        k1[4] = h*d2r(r,theta,d_r,d_theta, d_lambda, m, vol, g, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[5] = h*d2theta(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[6] = h*d2lambda(r,theta,d_r,d_theta,d_lambda,m,vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
        k1[7] = h*fluxo_massa(ro_h, d, P, P_h);
        k1[8] = h*dTh_dt(T, T_h, P, P_h, ro, d_r, r_b, r, ro_h, g, m, drb_dpm, alt_0, T0, m_b);

        //calculo k2
        k2[0] = h*drb_dt(r+k1[1]/2.0, d_r+k1[4]/2.0, T, P, r_b+k1[0]/2.0, ro_h, g, m+k1[7]/2.0, P_h, drb_dpm, alt_0, T0);
        k2[1] = h*(d_r+k1[4]/2.0);
        k2[2] = h*(d_theta+k1[5]/2.0);
        k2[3] = h*(d_lambda+k1[6]/2.0);
        k2[4] = h*d2r(r+k1[1]/2.0,theta+k1[2]/2.0,d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m+k1[7]/2.0,vol, g, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[5] = h*d2theta(r+k1[1]/2.0,theta+k1[2]/2.0,d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m+k1[7]/2.0,vol, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[6] = h*d2lambda(r+k1[1]/2.0, theta+k1[2]/2.0, d_r+k1[4]/2.0,d_theta+k1[5]/2.0, d_lambda+k1[6]/2.0,m+k1[7]/2.0,vol, ro, r_b+k1[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k2[7] = h*fluxo_massa(ro_h, d, P, P_h);
        k2[8] = h*dTh_dt(T, T_h+k1[8]/2.0, P, P_h, ro, d_r+k1[4]/2.0, r_b+k1[0]/2.0, r+k1[1]/2.0, ro_h, g, m, drb_dpm, alt_0, T0, m_b);

        //calculo de k3
        k3[0] = h*drb_dt(r+k2[1]/2.0, d_r+k2[4]/2.0, T, P, r_b+k2[0]/2.0, ro_h, g, m+k2[7]/2.0, P_h, drb_dpm, alt_0, T0);
        k3[1] = h*(d_r+k2[4]/2.0);
        k3[2] = h*(d_theta+k2[5]/2.0);
        k3[3] = h*(d_lambda+k2[6]/2.0);
        k3[4] = h*d2r(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m+k2[7]/2.0,vol, g, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[5] = h*d2theta(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m+k2[7]/2.0,vol, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[6] = h*d2lambda(r+k2[1]/2.0,theta+k2[2]/2.0,d_r+k2[4]/2.0,d_theta+k2[5]/2.0, d_lambda+k2[6]/2.0,m+k2[7]/2.0,vol, ro, r_b+k2[0]/2.0, vvento_r, vvento_theta, vvento_lambda);
        k3[7] = h*fluxo_massa(ro_h, d, P, P_h);
        k3[8] = h*dTh_dt(T, T_h+k2[8]/2.0, P, P_h, ro, d_r+k2[4]/2.0, r_b+k2[0]/2.0, r+k2[1]/2.0, ro_h, g, m, drb_dpm, alt_0, T0, m_b);

        //calculo de k4
        k4[0] = h*drb_dt(r+k3[1], d_r+k3[4], T, P, r_b+k3[0], ro_h, g, m+k3[7], P_h, drb_dpm, alt_0, T0);
        k4[1] = h*(d_r+k3[4]);
        k4[2] = h*(d_theta+k3[5]);
        k4[3] = h*(d_lambda+k3[6]);
        k4[4] = h*d2r(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m+k3[7],vol, g, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[5] = h*d2theta(r+k3[1],theta+k3[2],d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m+k3[7],vol, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[6] = h*d2lambda(r+k3[1],theta+k3[2], d_r+k3[4],d_theta+k3[5], d_lambda+k3[6],m+k3[7],vol, ro, r_b+k3[0], vvento_r, vvento_theta, vvento_lambda);
        k4[7] = h*fluxo_massa(ro_h, d, P, P_h);
        k4[8] = h*dTh_dt(T, T_h+k3[8], P, P_h, ro, d_r+k3[4], r_b+k3[0], r+k3[1], ro_h, g, m, drb_dpm, alt_0, T0, m_b);

        t += h;// O tempo � acrescido pelo passo h a cada loop
        cont_p++; //Contador de impress�o se soma

        //C�lculo das vari�veis
        r_b += (k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0])/6.0;
        r += (k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])/6.0;
        theta += (k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2])/6.0;
        lambda += (k1[3]+2.0*k2[3]+2.0*k3[3]+k4[3])/6.0;
        d_r += (k1[4]+2.0*k2[4]+2.0*k3[4]+k4[4])/6.0;
        d_theta += (k1[5]+2.0*k2[5]+2.0*k3[5]+k4[5])/6.0;
        d_lambda += (k1[6]+2.0*k2[6]+2.0*k3[6]+k4[6])/6.0;
        T_h += (k1[8]+2.0*k2[8]+2.0*k3[8]+k4[8])/6.0;

        //Vaz�o
        if ((alt>alt_flut-1000)&&(d_r>0)){
            m -= (k1[7]+2.0*k2[7]+2.0*k3[7]+k4[7])/6.0;
        }

        //empuxo = ro*vol*(g+(r*OMEGA*OMEGA*sin(theta)*sin(theta)));
        //peso = m*g;

        //Impress�o das vari�veis no arquivo
        if (cont_p==1000){ // A cada 1000 loops do while os valores das vari�veis s�o impressas no arquivo
            if(arq==NULL){ //Mensagem caso haja erro na abertura do arquivo
                printf("Problemas na abertura do arquivo\n");
                system("pause");
                exit(1);
            }
            printf("%13.8lf %13.8lf %13.8lf %13.8lf\n",t , r, d_r, m);//impress�o de dados desejados na tela
            //impress�o de dados no arquivo
            result = fprintf(arq,"%13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf %13.8lf\n",  t ,  r - R_T, d_r, m,  T, m_h, T_h);
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

