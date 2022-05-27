#include "prak6.h"

#define N 4		//���������� �����
#define K 3		//���������� �������� ���������
#define L 10	//���������� ���������� (�� �������� ��������) �����

int     main(void)
{
    double  a, b;					//���������� ������ �������
    size_t  M = (N - 1) * K + 1;	//����� ���������� �����
    double  *x;						//����������� ����� � ����� h
    double  h;						//��� ����������� ����� � ������ M �����

    a = -1.0;
    b = 1.0;
    h = (b - a) / (M - 1);

    //����������� ������ � ����� h
    x = x_gen(a, M, h);

    //���������� ������������� ��� ������� ���������� �����������
    slau(x, N, L, K);

    delete [] x;
    return 0;
}
