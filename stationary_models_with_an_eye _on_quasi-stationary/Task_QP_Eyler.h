#pragma once
#include "const.h"
#include "struct.h"

/// @brief Task_QP_Eiler - ����� ��� ������� ������ QP ��������� ������� ������
class Task_QP_Eyler
{
	// ���� ������
	// pipiline_parameters - c�������� ��������� ������������
	Pipiline_parameters pipiline_parameters;
	// oil_parameters - c�������� ��������� �����
	Oil_parameters oil_parameters;
	// tw - ����������� ���������� ������, ����������� ������ �������� ��� ������� �� �����
	double tw;
	// n - ��� - �� ����� ��������� �����
	int n;
	// h - ��� �� �� ���������� ��������� ����� [�]
	double h;
	// internal_diameter - ���������� ������� ����� [�]
	double internal_diameter;
	// pressure_previous - �������� �� ���������� �������� ��������� �������) [��]
	double pressure_previous;
	// pressure_current - �������� �� ������� ��������(������������ ��������) [��]

public:
	/// @brief ����������� ������
	/// @param pipiline_parameters - c�������� ��������� ������������
	/// @param oil_parameters - c�������� ��������� �����
	/// @param tw - ����������� ���������� ������, ����������� ������ �������� ��� ������� �� �����
	/// @param n - ��� - �� ����� ��������� �����
	/// @param h - ��� �� �� ���������� ��������� ����� [�]
	/// @param internal_diameter - ���������� ������� ����� [�]
	/// @param p_previous - �������� �� ���������� �������� ��������� �������) [��]
	Task_QP_Eyler(const Pipiline_parameters& pipiline_parameters,const Oil_parameters& oil_parameters,
		double tw, int n, double h, double internal_diameter, double p_previous);
	// ������ ������

	/// @brief ������ ������
	/// @param pipiline_parameters - c�������� ��������� ������������
	/// @param oil_parameters - c�������� ��������� �����
	/// @param tw - ����������� ���������� ������, ����������� ������ �������� ��� ������� �� �����
	/// @param n - ��� - �� ����� ��������� �����
	/// @param h - ��� �� �� ���������� ��������� ����� [�]
	/// @param internal_diameter - ���������� ������� ����� [�]
	/// @param p_previous - �������� �� ���������� �������� ��������� �������) [��]
	void setter(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
		double tw, int n, double h, double internal_diameter, double p_previous);
	
	/// @brief ������ ��� ������� ������ QP ��������� ������� ������
	/// @return pressure_current - �������� �� ������� ��������(������������ ��������)[��]
	double solver_eyler();

};       


 