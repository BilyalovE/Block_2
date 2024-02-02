#pragma once

#include "stationary_models_with_an_eye _on_quasi-stationary.cpp"

/// @brief Bernoulli_equation - ����� ��� ������� ����� �� ����� 2 - ���������� ������������ ������� 
/// � �������� �� �������������� (��������� ��������)
class Bernoulli_equation {
	// ���� ������
	// m_pipiline_parameters - ��������� ��������� ������������
	Pipiline_parameters m_pipiline_parameters;
	// m_oil_parameters - ��������� ��������� �����
	Oil_parameters m_oil_parameters;
	// m_hydraulic_resistance - ����������� ��������������_������������� (lambda)
	double m_hydraulic_resistance;
	// v - c������� ������� ����� [�/�]
	double m_v;
	// m_d - ���������� ������� �����[�]
	double m_d;

public:
	// ����������� ������ Bernoulli_equation
	Bernoulli_equation()

	/// @brief pressure_p0 - �����, �������������� �������� � ������ ������� ������������ [��]
	/// @param pipiline_parameters_XX - ��������� ��������� ������������
	/// @param oil_parameters_XX - ��������� ��������� �����
	/// @param lambda - ����������� ��������������� �������������
	/// @param v - �������� ������� ����� [�/�]
	/// @param d - ���������� ������� ����� [�]
	/// @return oil_parameters_XX.p0 - �������� � ������ ������� ������������ [��]
	double pressure_p0() {
		oil_parameters_XX.p0 = (oil_parameters_XX.ro * 9.81) * (oil_parameters_XX.pl * 1000000 / (oil_parameters_XX.ro * 9.81) - pipiline_parameters_XX.z0 + pipiline_parameters_XX.zl + (((lambda * pipiline_parameters_XX.l * 1000) / d) * pow(v, 2)) / (2 * 9.81)) * 0.000001;
		return oil_parameters_XX.p0;
	}
	/// @brief diameter - �����, �������������� ���������� ������� �����
	/// @param pipiline_parameters_XX - ��������� ��������� ������������
	/// @return d - ���������� ������� ����� [�]
	double diameter(Pipiline_parameters pipiline_parameters_XX) {
		double d = (pipiline_parameters_XX.D - 2 * pipiline_parameters_XX.ds) / 1000;
		return d;
	}

	/// @brief  relative_roughness - �����, �������������� ������������� ������������� �������������
	/// @param pipiline_parameters_XX - ��������� ��������� ������������
	/// @param d - ���������� ������� ����� [�]
	/// @return e o������������ ������������� �������������
	double relative_roughness(Pipiline_parameters pipiline_parameters_XX, double d) {
		double e = (pipiline_parameters_XX.delta / 1000) / d;
		return e;
	}

	/// @brief reynolds_number - �����, �������������� ����� ����������, ��� nu ���������� � ������� ��
	/// @param oil_parameters_XX - ��������� ��������� �����
	/// @param v - c������� ������� ����� � ������� ��
	/// @param d - ���������� ������� ����� � ������� ��
	/// @return Re - ����� ����������
	double reynolds_number(Oil_parameters oil_parameters_XX, double v, double d) {
		double Re = v * d / (oil_parameters_XX.nu * 0.000001);
		return Re;
	}

	/// @brief speed_flow - �����, �������������� �������� �� ��������� ������� �����
	/// @param oil_parameters_XX - ��������� ��������� �����
	/// @param d - ���������� ������� ����� [�]
	/// @return v - c������� ������� ����� � ������� ��
	double speed_flow(Oil_parameters oil_parameters_QP, double d) {
		double v = (4 * oil_parameters_QP.Q / 3600) / (3.14 * pow(d, 2));
		return v;
	}

	/// @brief speed_pressure - �����, �������������� �������� �� �������� � ������ PP
	/// @param pipiline_parameters_XX - ��������� ��������� ������������
	/// @param oil_parameters_XX - ��������� ��������� �����
	/// @param lambda - ����������� ��������������� �������������
	/// @param d - ���������� ������� ����� � ������� ��
	/// @param oil_parameters_XX.p0 - �������� � ������ ������� ������������ [���]
	/// @param oil_parameters_XX.pl - �������� � ����� ������� ������������ [���]
	/// @return v - c������� ������� ����� � ������� ��	
	double speed_pressure(Pipiline_parameters pipiline_parameters_XX, Oil_parameters oil_parameters_XX, double lambda, double d) {
		double v = pow((2 * 9.81 * d / pipiline_parameters_XX.l / 1000 * ((oil_parameters_XX.p0 - oil_parameters_XX.pl) * 1000000 / (oil_parameters_XX.ro * 9.81) + pipiline_parameters_XX.z0 - pipiline_parameters_XX.zl) / lambda), 0.5);
		return v;
	}

	/// @brief volume_flow - �����, �������������� �������� ������
	/// @param v - �������� ������� ����� � ������� ��
	/// @param d - ���������� ������� ����� [�]
	/// @return Q - �������� ������ [�^3/�]
	double volume_flow(double v, double d) {
		double Q = 3.14 * pow(d, 2) * v * 3600 / 4;
		return Q;
	}
};