
#include "biblioteca.h"
#include <math.h>

//Biblioteca com fun��es necess�rias para determina��o da trajet�ria de voo de um bal�o meteorol�gico com h�lio

double grav(double r)//C�lculo da acelera��o da gravidade real
{
    double C, g;
    //M_t = 5.972e24; //Massa da Terra
    //C_G = 6.67408e-11; //Constante gravitacional
    C = 3.9857605e14; // C = C_G*M_t
    g = C/(r*r);
return g;
}

// Determina��o de atmosfera padr�o internacional

double Temp(double r){//C�lculo da temperatura em fun��o da altitude em K
    double T0, grad_t, T, alt;
    alt = r-R_T; //Altitude em rela��o ao n�vel do mar
    if (alt <= 11000){ //C�lculo da temperatura na Troposfera
        T0 = 299.3734;
        grad_t = -4.9687e-3; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        T = T0 + grad_t*alt;
    }
    if ((11000 < alt) && (alt <= 25000)){ //C�lculo da temperatura na baixa estratosfera
        T0 = 244.7177;
        grad_t = -2.01555e-3;
        T = T0 + grad_t*(alt-11000);
    }
    if (alt > 25000){//C�lculo da temperatura na alta estratosfera
        T0 = 216.5;
        grad_t = 0.003; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        T = T0 + grad_t*(alt-25000);
    }
 return T;
}

double Pres(double r, double T, double g){ //C�lculo da press�o atmosf�rica local em Pa
    double T0, grad_t, P0, P, C;
    double alt = r - R_T; //Altitude em rela��o ao n�vel do mar
    if (alt <= 11000){ //C�lculo da press�o na troposfera
        T0 = 299.3734;
        P0 = 104036.51034;
        grad_t = -4.9687e-3; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
    }
    if ((11000 < alt) && (alt <= 25000)){ //C�lculo da press�o na baixa estratosfera
        T0 = 244.7177;
        P0 = 26039.4749;
        grad_t = -2.01555e-3; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
    }
    if (alt > 25000){ //C�lculo da press�o na alta estratosfera
        T0 = 216.5;
        P0 = 3268.6677;
        grad_t = 0.003; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        C = - g/(Ra*grad_t);
        P = P0 * pow((T/T0),C);
    }
 return P;
}

double dens(double r, double T, double P, double g){ // Densidade do ar por altitude (kg/m^3)
    double T0, grad_t, ro_0, C, P0, ro;
    double alt = r - R_T; //Altitude em rela��o ao n�vel do mar
    if (alt <= 11000){ //C�lculo da densidade na troposfera
        T0 = 299.3734;
        ro_0 = 1.225;
        grad_t = -4.9687e-3; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m (G3)
        C = -((g/(Ra*grad_t))+1);
        ro = ro_0 * pow((T/T0),C);
    }
    if ((11000 < alt) && (alt <= 25000)){ //C�lculo da densidade na baixa estratosfera
        T0 = 244.7177;
        P0 = 26037.064;
        ro_0 = P0/(Ra*T0);
        grad_t = -2.01555e-3; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        C = - g*(alt-11000)/(Ra*T0);
        ro = ro_0 * exp(C);
    }
    if (alt > 25000){//C�lculo da densidade na alta estratosfera
        T0 = 216.5;
        P0 = 2481;
        ro_0 = P0/(Ra*T0);
        grad_t = 0.003; //gradiente t�rmico espec�fico da camada atmosf�rica em K/m
        C = -((g/(Ra*grad_t))+1);
        ro = ro_0 * pow((T/T0),C);
    }
 return ro;
}

double dens_helio(double P, double T){ //Densidade local do h�lio
    double ro_h = P/(R_h*T);
    return ro_h;
}

double volume(double m_h, double ro_h){ //Volume local do bal�o
    double vol = m_h/ro_h;
    return vol;
}

//Determina��o do raio inicial do bal�o em fun��o da altitude
double raio_bi(double ro_h, double m_h)
{
    double r_bi = pow((3*m_h/(4*M_PI*ro_h)),1.0/3.0); //raio do bal�o conforme densidade de h�lio
    return r_bi;
}

//Determina��o do raio do bal�o pelo rela��o tens�o x deforma��o (1)
/*double raio_b(double P, double T, double m_h, double ro_h, double r_b, double vol, double e_0, double r_b0){
    //double r_bi = pow((3*m_h/(4*M_PI*ro_h)),1.0/3.0); //raio do bal�o conforme densidade inicial
    //double vol = 4.0*M_PI*pow(r_b,3.0)/3.0;
    double P_h = R_h*ro_h*T; //Press�o interna no bal�o
    //double drp = -3.205936e-6 - 4.886422e-5*(pow(0.999954, P_h));//derivada do raio do bal�o pela press�o interna
    double e = e_0*r_b0*r_b0/(r_b*r_b);
    double beta = 1.71583389E7 -4.81401394E8*pow(0.064193536,r_b); //r_b*r_b/(2.0*e*drp) + 3.0*r_b*r_b*r_b*(P_h-P)/(2.0*e_0*r_b0*r_b0); //M�dulo de young calculado (troposfera)
    double delta_v = (m_h/ro_h) - vol; //Varia��o do volume
    double r_bf = r_b*beta/(beta - (P_h-P)*r_b/2*e);// + (pow(3.0*delta_v/(4.0*M_PI),(1.0/3.0)));
    return r_bf;
}*/

//C�lculo de dr_b/dt para a determina��o do raio do bal�o no Runge-Kutta
double drb_dt(double r, double d_r, double T, double P, double r_b, double ro_h, double g, double m, double P_h, double drb_dpm)
{
    double alt, grad_t, drb_dtf;
    alt = r - R_T;
    if (alt <= 11000){
        grad_t = -4.9687e-3;
    }
    if ((alt > 11000) && (alt <= 25000)){
        grad_t = -2.01555e-3;
    }
    if (alt >25000){
        grad_t = 0.003;
    }
    drb_dtf = (drb_dpm*(d_r/T)*((P_h*grad_t) + (P*g/Ra)))/(1.0+(3.0*P_h*drb_dpm/r_b));
    return drb_dtf;
}

//C�lculo do arrasto gerado pelo bal�o na dire��o radial
double arrasto_r(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double v_modulo = sqrt(pow(v_r-vvento_r,2) + pow(v_theta-vvento_theta,2) + pow(v_lambda-vvento_lambda,2)); //M�dulo do vetor velocidade total
    double arrasto = -ro*(M_PI)*v_modulo*(v_r-vvento_r)*pow(r_b,2)*Ca/2;
return arrasto;
}

//C�lculo do arrasto gerado pelo bal�o na dire��o theta
double arrasto_theta(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double v_modulo = sqrt(pow(v_r-vvento_r,2) + pow(v_theta-vvento_theta,2) + pow(v_lambda-vvento_lambda,2)); //M�dulo do vetor velocidade total
    double arrasto = -ro*(M_PI)*v_modulo*(v_theta-vvento_theta)*pow(r_b,2)*Ca/2;
return arrasto;
}

//C�lculo do arrasto gerado pelo bal�o na dire��o lambda
double arrasto_lambda(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double v_modulo = sqrt(pow(v_r-vvento_r,2) + pow(v_theta-vvento_theta,2) + pow(v_lambda-vvento_lambda,2)); //M�dulo do vetor velocidade total
    double arrasto = -ro*(M_PI)*v_modulo*(v_lambda-vvento_lambda)*pow(r_b,2)*Ca/2;
return arrasto;
}

//Acelera��o na dire��o radial
double d2r(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double g, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double s = sin(theta);
    double arrasto = arrasto_r(d_r, d_theta*r, d_lambda*r*s, vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
    double d2r_r = (ro*vol*(g - OMEGA*OMEGA*r*s*s))/m - g + arrasto/m + OMEGA*r*s*s*(2*d_lambda + OMEGA) + r*d_theta*d_theta + r*s*s*d_lambda*d_lambda;
    return d2r_r;
}

//Acelera��o na dire��o theta
double d2theta(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)
{
    double s = sin(theta);
    double arrasto = arrasto_theta(d_r, d_theta*r, d_lambda*r*s, vol, ro, r_b, vvento_r, vvento_theta, vvento_lambda);
    double sc = s*cos(theta);
    double d2theta_r = (-ro*vol*OMEGA*OMEGA*sc)/m + arrasto/(r*m) + OMEGA*sc*(2*d_lambda + OMEGA) - 2*d_r*d_theta/r + sc*d_lambda*d_lambda;
    return d2theta_r;
}

//Acelera��o na dire��o lambda
double d2lambda(double r, double theta, double d_r,double d_theta, double d_lambda,double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda)//Acelera��o na dire��o lambda
{
    double s = sin(theta);
    double c = cos(theta);
    double arrasto = arrasto_lambda(d_r, d_theta*r, d_lambda*r*s, vol, ro, r_b,vvento_r, vvento_theta, vvento_lambda);
    //double dvlambda_r = arrasto/(r*s*m) - 2*OMEGA*(v_r/r + v_theta*c/s) - 2*v_r*v_lambda/r - 2*v_theta*v_lambda*c/s;
    double d2lambda_r = arrasto/(r*s*m) - 2*d_r/r*(OMEGA + d_lambda) - 2*d_theta*c/s*(OMEGA + d_lambda);
    return d2lambda_r;
}

