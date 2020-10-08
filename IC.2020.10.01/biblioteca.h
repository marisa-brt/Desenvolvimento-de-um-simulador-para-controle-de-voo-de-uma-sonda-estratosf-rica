

//Biblioteca com fun��es necess�rias para determina��o da trajet�ria de voo de um bal�o meteorol�gico com h�lio

#define R_T 6371000 //Raio da Terra
#define OMEGA 7.292115900231275e-5 // Velocidade angular da Terra
#define R_h 2077.1 //Constante individual do h�lio
#define Ca 0.34037379532479679571 //Constante de arrasto
#define Ra 287.54 // Constante universal dos gases tomando a massa molar do ar constante


double grav(double r); //Fun��o acelera��o da gravidade local em m/s^2
double update_atm(double alt_0, double T0, double r, double r_ant, double g, double *T, double *p, double *rho);//Determina as propriedades termodin�micas da atmosfera na posicao radial r
double Temp(double r); //Fun��o temperatura por altitude em K
double Pres(double r, double T, double g); //Fun��o press�o por altitude em Pa
double dens(double r, double T, double P, double g); //Fun��o densidade por altitude em kg/m^3
double dens_helio(double P, double T);//Densidade local do h�lio
double volume(double m_h, double ro_h); //Volume local do bal�o
double raio_bi(double ro_h, double m_h); //Fun��o raio do bal�o inicial
double drb_dt(double r, double d_r, double T, double P, double r_b, double ro_h, double g, double m, double P_h, double drb_dpm);
double arrasto_r(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Fun��o arrasto gerado pelo bal�o em r
double arrasto_theta(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Fun��o arrasto gerado pelo bal�o em theta
double arrasto_lambda(double v_r,double v_theta, double v_lambda, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Fun��o arrasto gerado pelo bal�o em lambda
double d2r(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double g, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Acelera��o radial
double d2theta(double r, double theta, double d_r, double d_theta, double d_lambda, double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Acelera��o em theta
double d2lambda(double r, double theta, double d_r,double d_theta, double d_lambda,double m, double vol, double ro, double r_b, double vvento_r, double vvento_theta, double vvento_lambda); //Acelera��o em lambda
