#include "biblioteca.h"
#include <math.h>

//Biblioteca com funções necessárias para determinação da trajetória de voo de um balão meteorológico com hélio

double grav(double r)//Cálculo da aceleração da gravidade real
{
    double C, g;
    //M_t = 5.972e24; //Massa da Terra
    //C_G = 6.67408e-11; //Constante gravitacional
    C = 3.9857605e14; // C = C_G*M_t
    g = C/(r*r);
return g;
}

double update_atm(double alt_0, double T0, double r, double r_ant, double g, double *T, double *p, double *rho){

  //Determina as propriedades termodinâmicas da atmosfera na posicao radial r
  //T em K
  //p em Pa
  //rho em kg/m^3

  //alt_0 = altitude do local do lancamento
  //T0 = temperatura no local do lancamento
  //r = coordenada radial
  //r_ant = coordenada radial na iteração anterior
  //g = aceleracao da gravidade na coordenada r

  double grad_t1 = -5.830194e-3; //gradiente térmico na Troposfera (K/m)
  double grad_t3 = 0.003; //gradiente térmico na alta estratosfera (K/m)
  double T2 = 216.5; //temperatura da 2a camada

  double alt_12 = alt_0 + (T2-T0)/grad_t1;
  //altitude de transicao entre a 1a e a 2a camadas

  double alt = r-R_T; //Altitude em relação ao nível do mar

  double T_ant = *T;
  double p_ant = *p;

  double pot1 = -g/Ra/grad_t1;
  double pot2 = -g/Ra/T2;
  double pot3 = -g/Ra/grad_t3;

  if(alt <= alt_12){
    //Cálculo da temperatura na Troposfera
    *T = T0 + grad_t1*(alt - alt_0);
    *p = p_ant*pow((*T)/T_ant,pot1);
    *rho = *p/Ra/(*T);
  }
  else if((alt_12 <= alt) && (alt <= 25000)){
    //Cálculo da temperatura na baixa estratosfera
    *T = T2;
    *p = p_ant*exp(pot2*(r-r_ant));
    *rho = *p/Ra/(*T);
  }
  else{
    //Cálculo da temperatura na alta estratosfera
    *T = T2 + grad_t3 * (alt - 25000);
    *p = p_ant * pow((*T)/T_ant,pot3);
    *rho = *p /Ra/(*T);
  }
  //return T;
}


// Determinação de atmosfera padrão internacional

double Temp(double r){//Cálculo da temperatura em função da altitude em K
    double T0, grad_t, T, alt;
    alt = r-R_T; //Altitude em relação ao nível do mar
    if (alt <= 11000){ //Cálculo da temperatura na Troposfera
        T0 = 299.3734;
        grad_t = -4.9687e-3; //gradiente térmico específico da camada atmosférica em K/m
        T = T0 + grad_t*alt;
    }
    if ((11000 < alt) && (alt <= 25000)){ //Cálculo da temperatura na baixa estratosfera
        T0 = 244.7177;
        grad_t = -2.01555e-3;
        T = T0 + grad_t*(alt-11000);
    }
    if (alt > 25000){//Cálculo da temperatura na alta estratosfera
        T0 = 216.5;
        grad_t = 0.003; //gradiente térmico específico da camada atmosférica em K/m
        T = T0 + grad_t*(alt-25000);
    }
 return T;
}

double Pres(double r, double T, double g){ //Cálculo da pressão atmosférica local em Pa
    double T0, grad_t, P0, P, C;
    double alt = r - R_T; //Altitude em relação ao nível do mar
    if (alt <= 11000){ //Cálculo da pressão na troposfera
        T0 = 299.3734;
        P0 = 104036.51034;
        grad_t = -4.9687e-3; //gradiente térmico específico da camada atmosférica em K/m
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
    }
    if ((11000 < alt) && (alt <= 25000)){ //Cálculo da pressão na baixa estratosfera
        T0 = 244.7177;
        P0 = 26039.4749;
        grad_t = -2.01555e-3; //gradiente térmico específico da camada atmosférica em K/m
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
    }
    if (alt > 25000){ //Cálculo da pressão na alta estratosfera
        T0 = 216.5;
        P0 = 3268.6677;//
        grad_t = 0.003; //gradiente térmico específico da camada atmosférica em K/m
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
    }
 return P;
}

double dens(double r, double T, double P, double g){ // Densidade do ar por altitude (kg/m^3)
    double T0, grad_t, ro_0, C, P0, ro;
    double alt = r - R_T; //Altitude em relação ao nível do mar
    if (alt <= 11000){ //Cálculo da densidade na troposfera
        T0 = 299.3734;
        ro_0 = 1.225;
        grad_t = -4.9687e-3; //gradiente térmico específico da camada atmosférica em K/m (G3)
        C = -((g/(Ra*grad_t))+1);
        ro = ro_0 * pow((T/T0),C);
    }
    if ((11000 < alt) && (alt <= 25000)){ //Cálculo da densidade na baixa estratosfera
        T0 = 244.7177;
        P0 = 26039.4749;
        ro_0 = P0/(Ra*T0);
        grad_t = -2.01555e-3; //gradiente térmico específico da camada atmosférica em K/m
        C = -((g/(Ra*grad_t))+1);
        ro = ro_0 * pow((T/T0),C);
    }
    if (alt > 25000){//Cálculo da densidade na alta estratosfera
        T0 = 216.5;
        P0 = 3268.6677;
        ro_0 = P0/(Ra*T0);
        grad_t = 0.003; //gradiente térmico específico da camada atmosférica em K/m
        C = -((g/(Ra*grad_t))+1);
        ro = ro_0 * pow((T/T0),C);
    }
 return ro;
}

double dens_helio(double P, double T){ //Densidade local do hélio
    double ro_h = P/(R_h*T);
    return ro_h;
}

double volume(double m_h, double ro_h){ //Volume local do balão
    double vol = m_h/ro_h;
    return vol;
}

//Determinação do raio inicial do balão em função da altitude
double raio_bi(double ro_h, double m_h)
{
    double r_bi = pow((3*m_h/(4*M_PI*ro_h)),1.0/3.0); //raio do balão conforme densidade de hélio
    return r_bi;
}

//Cálculo de dr_b/dt para a determinação do raio do balão no Runge-Kutta
double drb_dt(double r, double d_r, double T, double P, double r_b, double ro_h, double g, double m, double P_h, double drb_dpm, double alt_0, double T0)
{
    double alt, grad_t, drb_dtf;
    double grad_t1 = -5.830194e-3; //gradiente térmico na Troposfera (K/m)
    double grad_t3 = 0.003; //gradiente térmico na alta estratosfera (K/m)
    double T2 = 216.5; //temperatura da 2a camada

    double alt_12 = alt_0 + (T2-T0)/grad_t1;
    //altitude de transicao entre a 1a e a 2a camadas

    alt = r - R_T;
    if (alt <= alt_12){
        grad_t = grad_t1;
    }
    if ((alt > alt_12) && (alt <= 25000)){
        grad_t = 0;
    }
    if (alt >25000){
        grad_t = grad_t3;
    }
    drb_dtf = (drb_dpm*(d_r/T)*((P_h*grad_t) + (P*g/Ra)))/(1.0+(3.0*P_h*drb_dpm/r_b));
    return drb_dtf;
}

//Cálculo do arrasto gerado pelo balão na direção radial
double arrasto_r(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double v_modulo = sqrt(pow(v_r-vvento_r,2) + pow(v_theta-vvento_theta,2) + pow(v_lambda-vvento_lambda,2)); //Módulo do vetor velocidade total
    double arrasto = -ro*(M_PI)*v_modulo*(v_r-vvento_r)*pow(r_b,2)*Ca/2;
return arrasto;
}

//Cálculo do arrasto gerado pelo balão na direção theta
double arrasto_theta(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double v_modulo = sqrt(pow(v_r-vvento_r,2) + pow(v_theta-vvento_theta,2) + pow(v_lambda-vvento_lambda,2)); //Módulo do vetor velocidade total
    double arrasto = -ro*(M_PI)*v_modulo*(v_theta-vvento_theta)*pow(r_b,2)*Ca/2;
return arrasto;
}

//Cálculo do arrasto gerado pelo balão na direção lambda
double arrasto_lambda(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double v_modulo = sqrt(pow(v_r-vvento_r,2) + pow(v_theta-vvento_theta,2) + pow(v_lambda-vvento_lambda,2)); //Módulo do vetor velocidade total
    double arrasto = -ro*(M_PI)*v_modulo*(v_lambda-vvento_lambda)*pow(r_b,2)*Ca/2;
return arrasto;
}

//Aceleração na direção radial
double d2r(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double g, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double s = sin(theta);
    double arrasto = arrasto_r(d_r, d_theta*r, d_lambda*r*s, vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
    double d2r_r = (ro*vol*(g - OMEGA*OMEGA*r*s*s))/m - g + arrasto/m + OMEGA*r*s*s*(2*d_lambda + OMEGA) + r*d_theta*d_theta + r*s*s*d_lambda*d_lambda;
    return d2r_r;
}

//Aceleração na direção theta
double d2theta(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double s = sin(theta);
    double arrasto = arrasto_theta(d_r, d_theta*r, d_lambda*r*s, vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
    double sc = s*cos(theta);
    double d2theta_r = (-ro*vol*OMEGA*OMEGA*sc)/m + arrasto/(r*m) + OMEGA*sc*(2*d_lambda + OMEGA) - 2*d_r*d_theta/r + sc*d_lambda*d_lambda;
    return d2theta_r;
}

//Aceleração na direção lambda
double d2lambda(double r, double theta, double d_r,double d_theta, double d_lambda,double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)//Aceleração na direção lambda
{
    double s = sin(theta);
    double c = cos(theta);
    double arrasto = arrasto_lambda(d_r, d_theta*r, d_lambda*r*s, vol, ro, r_b,vvento_r, vvento_theta, vvento_lambda);
    //double dvlambda_r = arrasto/(r*s*m) - 2*OMEGA*(v_r/r + v_theta*c/s) - 2*v_r*v_lambda/r - 2*v_theta*v_lambda*c/s;
    double d2lambda_r = arrasto/(r*s*m) - 2*d_r/r*(OMEGA + d_lambda) - 2*d_theta*c/s*(OMEGA + d_lambda);
    return d2lambda_r;
}

//Fluxo de massa na válvula
double fluxo_massa(double ro_h, double d, double P, double P_h)
{
    double dm_ds = M_PI*d*d*0.25*sqrt(2.0*(P_h-P)/(K*ro_h));
    return dm_ds;
}

