#pragma once
#include "struct.h"
#include "Task_QP_Eyler.h"
#include <iomanip>
#include <fixed/fixed.h>
class Task_PP_Newton_Eyler : public Task_QP_Eyler, public fixed_system_t<1>
{
	Oil_parameters m_oil_parameters;
	Pipeline_parameters m_pipeline_parameters;
	double v;

	using fixed_system_t<1>::var_type;
public:

	// Конструктор класса
	Task_PP_Newton_Eyler(Pipeline_parameters pipeline_parameters, Oil_parameters oil_parameters, double v);
	
	
	/// @brief residuals - ������� �������
	/// @param - v - ������� �������� (��������, [�/�])
	var_type residuals(const var_type& x) {
		double count_pressure_p0 = solver_eyler(x);
		return (m_oil_parameters.p0 - count_pressure_p0);
	}

	double solver_newton_rafson();
};


