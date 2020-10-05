//Biblioteca com funções necessárias para determinação da trajetória de voo de um balão meteorológico com hélio

#define R_T 6371000 //Raio da Terra
#define OMEGA 7.292115900231275e-5 // Velocidade angular da Terra
#define R_h 2077.1 //Constante individual do hélio
#define Ca 0.34037379532479679571 //Constante de arrasto
#define Ra 287.54 // Constante universal dos gases tomando a massa molar do ar constante
//#define Elastic 0.1e6 //Módulo de Young médio para o látex
//#define e 0.001 //Espessura


double grav(double r); //Função aceleração da gravidade local em m/s^2
double Temp(double r); //Função temperatura por altitude em K
double Pres(double r, double T, double g); //Função pressão por altitude em Pa
double dens(double r, double T, double P, double g); //Função densidade por altitude em kg/m^3
double update_atm(double alt_0, double T0, double r, double r_ant, double g,
		  double *T, double *p, double *rho);
double dens_helio(double P, double T);//Densidade local do hélio
double volume(double m_h, double ro_h); //Volume local do balão
double raio_bi(double ro_h, double m_h); //Função raio do balão inicial
double drb_dt(double r, double d_r, double T, double P, double r_b, double ro_h, double g, double m, double P_h, double drb_dpm);
//double raio_b(double P, double T, double m_h, double ro_h, double r_b, double vol, double e_0, double r_b0);//Função parao raio do balão em função da altitude e da deformação
double arrasto_r(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Função arrasto gerado pelo balão em r
double arrasto_theta(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Função arrasto gerado pelo balão em theta
double arrasto_lambda(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Função arrasto gerado pelo balão em lambda
double d2r(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double g, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Aceleração radial
double d2theta(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Aceleração em theta
double d2lambda(double r, double theta, double d_r,double d_theta, double d_lambda,double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Aceleração em lambda

