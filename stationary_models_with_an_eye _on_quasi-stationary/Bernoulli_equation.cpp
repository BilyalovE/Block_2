
#include "Bernoulli_equation.h"

/// @brief Bernoulli_equation - конструктор класса по умолчанию
/// @param pipeline_parameters - структура параметров трубопровода
/// @param oil_parameters - структура параметров нефти
Bernoulli_equation::Bernoulli_equation(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters, double hydraulic_resistance)
{
	setter(pipeline_parameters, oil_parameters, hydraulic_resistance);

}



Bernoulli_equation::Bernoulli_equation(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters)
{
	m_pipeline_parameters = pipeline_parameters;
	m_oil_parameters = oil_parameters;
	m_d = m_pipeline_parameters.get_inner_diameter();
}

// Методы класса
void Bernoulli_equation::setter(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters, double hydraulic_resistance)
{
	m_pipeline_parameters = pipeline_parameters;
	m_oil_parameters = oil_parameters;
	m_hydraulic_resistance = hydraulic_resistance;
	m_d = m_pipeline_parameters.get_inner_diameter();
	m_Re = m_pipeline_parameters.reynolds_number();
	m_relative_roughness = m_pipeline_parameters.get_relative_roughness();
	m_v = m_pipeline_parameters.speed_flow();
}

double Bernoulli_equation::pressure_p0() {

	m_p0 = (m_oil_parameters.ro * k_g) * (m_oil_parameters.pl / (m_oil_parameters.ro * k_g)
		- m_pipeline_parameters.z0 + m_pipeline_parameters.zl
		+ (((m_hydraulic_resistance * m_pipeline_parameters.l) / m_d) * pow(m_v, 2)) / (2 * k_g));
	return m_p0;
}

double Bernoulli_equation::internal_diameter() {
	return m_d = (m_pipeline_parameters.D - 2 * m_pipeline_parameters.ds);
}

double Bernoulli_equation::relative_roughness() {
	return m_relative_roughness = (m_pipeline_parameters.delta) / m_d;
}

double Bernoulli_equation::reynolds_number() {
	return m_Re = m_v * m_d / (m_pipeline_parameters.nu);
}

double Bernoulli_equation::reynolds_number(double m_v) {
	return m_Re = m_v * m_d / (m_pipeline_parameters.nu);
}

double Bernoulli_equation::speed_flow() {
	return m_v = (4 * m_pipeline_parameters.Q) / (k_pi * pow(m_d, 2));
}

double Bernoulli_equation::speed_pressure() {
	m_v = pow((2 * k_g * m_d / m_pipeline_parameters.l * ((m_oil_parameters.p0 - m_oil_parameters.pl) / (m_oil_parameters.ro * k_g) + m_pipeline_parameters.z0 - m_pipeline_parameters.zl) / m_hydraulic_resistance), 0.5);
	return m_v;
}



double Bernoulli_equation::speed_pressure(double m_hydraulic_resistance) {
	m_v = pow((2 * k_g * m_d / m_pipeline_parameters.l * ((m_oil_parameters.p0 - m_oil_parameters.pl) / (m_oil_parameters.ro * k_g) + m_pipeline_parameters.z0 - m_pipeline_parameters.zl) / m_hydraulic_resistance), 0.5);
	return m_v;
}

double Bernoulli_equation::volume_flow() {
	return m_Q = k_pi * pow(m_d, 2) * m_v / 4;
}

double Bernoulli_equation::volume_flow(double m_v) {
	return m_Q = k_pi * pow(m_d, 2) * m_v / 4;
}