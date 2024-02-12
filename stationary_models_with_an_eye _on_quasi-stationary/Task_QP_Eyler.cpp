#include "Task_QP_Eyler.h"

Task_QP_Eyler::Task_QP_Eyler(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
	double tw, int n, double h, double internal_diameter, double p_previous)
{
	setter(pipiline_parameters, oil_parameters, tw, n, h, internal_diameter, p_previous);
}

// Методы класса

void Task_QP_Eyler::setter(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
	double tw, int n, double h, double internal_diameter, double pressure_previous) {
	this->pipiline_parameters = pipiline_parameters;
	this->oil_parameters = oil_parameters;
	this->tw = tw;
	this->n = n;
	this->h = h;
	this->internal_diameter = internal_diameter;
	this->pressure_previous = pressure_previous;
}
double Task_QP_Eyler::solver_eyler()
{
	// pressure_current - давление на текущей итерации(рассчитанное значение) [Па]
	double pressure_current;
	for (size_t i = 1; i <= n; i++) {
		pressure_current = (pressure_previous - h * (-4 * tw / internal_diameter - oil_parameters.ro * k_g * (pipiline_parameters.zl - pipiline_parameters.z0) / (n - 1) / h));
		pressure_previous = pressure_current;
	}
	return pressure_current;
}

