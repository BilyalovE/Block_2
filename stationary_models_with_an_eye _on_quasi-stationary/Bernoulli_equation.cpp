
#include "Bernoulli_equation.h"

/// @brief Bernoulli_equation - конструктор класса по умолчанию
/// @param pipiline_parameters - cтруктура парметров трубопровода
/// @param oil_parameters - структура парметров нефти
Bernoulli_equation::Bernoulli_equation(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters)
{
	setter1(pipiline_parameters, oil_parameters);

}

/// @brief Bernoulli_equation - конструктор класса
/// @param pipiline_parameters - труктура парметров трубопровода
/// @param oil_parameters - структура парметров нефти
/// @param hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
/// @param v - cкорость течения нефти [м/с]
/// @param d - внутренний диаметр трубы [м]
Bernoulli_equation::Bernoulli_equation(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
	double hydraulic_resistance, double v, double d)

{
	setter2(pipiline_parameters, oil_parameters, hydraulic_resistance, v, d);

}

// Методы класса
void Bernoulli_equation::setter1(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters)
{
	m_pipiline_parameters = pipiline_parameters;
	m_oil_parameters = oil_parameters;
}

void Bernoulli_equation::setter2(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
	double& hydraulic_resistance, const double& v, double& d)
{
	m_pipiline_parameters = pipiline_parameters;
	m_oil_parameters = oil_parameters;
	m_hydraulic_resistance = hydraulic_resistance;
	m_v = v;
	m_d = d;
}

double Bernoulli_equation::pressure_p0() {
	m_p0 = (m_oil_parameters.ro * k_g) * (m_oil_parameters.pl / (m_oil_parameters.ro * k_g)
		- m_pipiline_parameters.z0 + m_pipiline_parameters.zl
		+ (((m_hydraulic_resistance * m_pipiline_parameters.l) / m_d) * pow(m_v, 2)) / (2 * k_g));
	return m_p0;
}

double Bernoulli_equation::internal_diameter() {
	return m_d = (m_pipiline_parameters.D - 2 * m_pipiline_parameters.ds);
}

double Bernoulli_equation::relative_roughness() {
	return m_relative_roughness = (m_pipiline_parameters.delta) / m_d;
}

double Bernoulli_equation::reynolds_number() {
	return m_Re = m_v * m_d / (m_oil_parameters.nu);
}

double Bernoulli_equation::reynolds_number(double m_v) {
	return m_Re = m_v * m_d / (m_oil_parameters.nu);
}

double Bernoulli_equation::speed_flow() {
	return m_v = (4 * m_oil_parameters.Q) / (k_pi * pow(m_d, 2));
}

double Bernoulli_equation::speed_pressure() {
	m_v = pow((2 * k_g * m_d / m_pipiline_parameters.l * ((m_oil_parameters.p0 - m_oil_parameters.pl) / (m_oil_parameters.ro * k_g) + m_pipiline_parameters.z0 - m_pipiline_parameters.zl) / m_hydraulic_resistance), 0.5);
	return m_v;
}

double Bernoulli_equation::speed_pressure(double m_hydraulic_resistance) {
	m_v = pow((2 * k_g * m_d / m_pipiline_parameters.l * ((m_oil_parameters.p0 - m_oil_parameters.pl) / (m_oil_parameters.ro * k_g) + m_pipiline_parameters.z0 - m_pipiline_parameters.zl) / m_hydraulic_resistance), 0.5);
	return m_v;
}

double Bernoulli_equation::volume_flow() {
	return m_Q = k_pi * pow(m_d, 2) * m_v / 4;
}

double Bernoulli_equation::volume_flow(double m_v) {
	return m_Q = k_pi * pow(m_d, 2) * m_v / 4;
}