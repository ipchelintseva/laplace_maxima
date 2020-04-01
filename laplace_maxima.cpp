// ********************************************************
// * (c) 2020 Irina Pchelintseva, Alexander Pchelintsev,  *
// *          Yuriy Litovka, irina_yu_10@mail.ru          *
// * Решение уравнения Лапласа для поиска распределения   *
// * потенциала в гальванической ванной и расчёта толщины *
// * покрытия катода (никелирование)                      *
// ********************************************************

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

#include <cstdio>
#include <cstdlib>

#define Nx      15
#define Ny      15
#define hx      0.3
#define hy      0.3
#define Ua      5
#define kapa    0.515
#define E       1.09
#define rho     8.902
#define delta_t 0.5

#define min_i_a Nx/4 + 1
#define max_i_a 3*Nx/4
#define min_i_k Nx/3 + 1
#define max_i_k 2*Nx/3

#define FORMAT  "png"

void GetStringValue(ifstream &f, string &res)
{
	res = "";
	char c;
	while(f.get(c))
		if(c == '=')
			break;
	if(c == EOF)
		return;
	while(f.get(c))
	{
		if(c == ',' || c == ']')
			break;
		if(c != '\n' && c != ' ')
			res += c;
	}
}

void create_laplas_df(ofstream &newton)
{
	for(int i = 2; i <= Nx - 1; i++)
		for(int j = 2; j <= Ny - 1; j++)
		{
			newton << "(phi_" << i - 1 << "_" << j << "-2*phi_" << i << "_" << j <<
			          "+phi_" << i+1 << "_" << j <<
			          ")/" << hx*hx;
			newton << "+(phi_" << i << "_" << j-1 << "-2*phi_" << i << "_" << j <<
			          "+phi_" << i << "_" << j+1 << ")/" << hy*hy << "," << endl;
		}
}

void create_gran_izolyator_hor(ofstream &newton)
{
	for(int j = 2; j <= Ny - 1; j++)
	{
		newton << "phi_2_" << j << "-" << "phi_1_" << j << "," << endl;
		newton << "phi_" << Nx << "_" << j << "-" << "phi_" << Nx-1 <<
		          "_" << j << "," << endl;
	}
}

void create_gran_izolyator_vert_left1(ofstream &newton)
{
	for(int i = 1; i <= min_i_a - 1; i++)
		newton << "phi_" << i << "_2-" << "phi_" << i << "_1," << endl;
}

void create_gran_F1(ofstream &newton)
{
	for(int i = min_i_a; i <= max_i_a; i++)
		newton << "phi_" << i << "_1+Fa(" << kapa << "*(" <<
		          "phi_" << i << "_1-" << "phi_" << i << "_2)/" <<
		          hy << ")-" << Ua << "," << endl;
}

void create_gran_izolyator_vert_left2(ofstream &newton)
{
	for(int i = max_i_a + 1; i <= Nx; i++)
		newton << "phi_" << i << "_2-" << "phi_" << i << "_1," << endl;
}

void create_gran_izolyator_vert_right1(ofstream &newton)
{
	for(int i = 1; i <= min_i_k - 1; i++)
		newton << "phi_" << i << "_" << Ny << "-" <<
		          "phi_" << i << "_" << Ny-1 << "," << endl;
}

void create_gran_F2(ofstream &newton)
{
	for(int i = min_i_k; i <= max_i_k; i++)
		newton << "phi_" << i << "_" << Ny << "+" << "Fk(" <<
		          kapa << "*(phi_" << i << "_" << Ny-1 << "-" <<
		          "phi_" << i << "_" << Ny << ")/" << hy << ")," << endl;
}

void create_gran_izolyator_vert_right2(ofstream &newton)
{
	for(int i = max_i_k + 1; i <= Nx; i++)
	{
		newton << "phi_" << i << "_" << Ny << "-" << "phi_" << i << "_" <<
		          Ny-1;
		if(i != Nx)
			newton << "," << endl;
		else
			newton << endl;
	}
}

void create_vars(ofstream &newton)
{
	for(int i = 1; i <= Nx; i++)
	{
		for(int j = 1; j <= Ny; j++)
		{
			newton << "phi_" << i << "_" << j;
			if(i != Nx || j != Ny)
				newton << ",";
		}
		newton << endl;
	}
	newton << "]," << endl << "[" << endl;
}

void create_init_vars(ofstream &newton)
{
	for(int i = 1; i <= Nx; i++)
	{
		for(int j = 1; j <= Ny; j++)
		{
			newton << Ua * (j - Ny) / (double)(1 - Ny);
			if(i != Nx || j != Ny)
				newton << ",";
		}
		newton << endl;
	}
	newton << "]" << endl;
}

int main()
{
	ofstream newton("newton.wxm");
	newton << "/* [wxMaxima batch file version 1] [ DO NOT EDIT BY HAND! ]*/" <<
	          endl;
	newton << "/* [wxMaxima: input   start ] */" << endl;
	newton << "Fa(ia):=-4.267*ia^2+5.867*ia$" << endl;
	newton << "Fk(ik):=0.883*ik^2-2.242*ik$" << endl;
	newton << "display2d:false$" << endl;
	newton << "load(mnewton);" << endl;
	newton << "newtonepsilon:10^-2$" << endl;
	newton << "mnewton(" << endl << "[" << endl;
	create_laplas_df(newton);
	create_gran_izolyator_hor(newton);
	create_gran_izolyator_vert_left1(newton);
	create_gran_F1(newton);
	create_gran_izolyator_vert_left2(newton);
	create_gran_izolyator_vert_right1(newton);
	create_gran_F2(newton);
	create_gran_izolyator_vert_right2(newton);
	newton << "]," << endl << "[" << endl;
	create_vars(newton);
	create_init_vars(newton);
	newton << ");" << endl;
	newton << "/* [wxMaxima: input   end   ] */" << endl;
	newton.close();

	system("maxima < newton.wxm > res_maxima.txt");

	double **phi = new double* [Nx + 1];
	for(int i = 1; i <= Nx; i++)
		phi[i] = new double [Ny + 1];
	ifstream outm("res_maxima.txt");

	string s;
	ofstream phi_file("phi.txt");
	for(int i = 1; i <= Nx; i++)
	{
		for(int j = 1; j <= Ny; j++)
		{
			GetStringValue(outm, s);
			phi[i][j] = atof(s.c_str());
			phi_file << setiosflags(ios::right) << setw(8) << setprecision(3) <<
			            setfill(' ') << phi[i][j];
		}
		phi_file << endl;
	}
	outm.close();
	phi_file.close();

	double *delta = new double [Nx + 1];
	for(int i = 1; i <= Nx; i++)
		if(i < min_i_k || i > max_i_k)
			delta[i] = 0;
		else
			delta[i] = E * (0.1 * kapa) * (0.1 *(phi[i][Ny - 1] - phi[i][Ny])/(hy*rho)) * delta_t;
	for(int i = 1; i <= Nx; i++)
		delete []phi[i];
	delete []phi;

	ofstream graph_delta("graphic_delta.txt");
	graph_delta << "set term " << FORMAT << " size 1100, 700 font 20" <<
	               "\nset output \""<< "graphic_delta." << FORMAT << "\"\n";
	graph_delta << "set xrange [" << (min_i_k - 1) * hx << ":" <<
	               (max_i_k - 1) * hx << "]\n";
	graph_delta << "set xlabel \"x\"\n";
	graph_delta << "set ylabel \"delta\"\n";
	graph_delta << "plot \"-\" with points pt 7 ps 2 lc 3 title \"\",\\\n";
	graph_delta << "     \"-\" with lines lw 2 lc 1 title \"\" smooth csplines\n";
	for(int i = min_i_k; i <= max_i_k; i++)
			graph_delta << (i - 1) * hx << " " << 10000 * delta[i] <<  endl;
	graph_delta << "e\n";
	for(int i = min_i_k; i <= max_i_k; i++)
			graph_delta << (i - 1) * hx << " " << 10000 * delta[i] <<  endl;
	graph_delta.close();

	double delta_min = delta[min_i_k];
	for(int i = min_i_k + 1; i <= max_i_k; i++)
		if(delta[i] < delta_min)
			delta_min = delta[i];
	double R = 0;
	for(int i = min_i_k; i <= max_i_k; i++)
		R += delta[i] - delta_min;
	R /= delta_min * (max_i_k - min_i_k + 1);
	cout << "R = " << R << endl;

	delete []delta;
	return 0;
}
