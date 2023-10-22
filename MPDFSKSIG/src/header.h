#pragma once

#include <iostream>
#include <vector>
#include <time.h>
#include <string>

#include <thread>
#include "oneapi/tbb.h"
//#include <omp_llvm.h>
#include <omp.h>

#include <boost/multiprecision/cpp_int.hpp>

extern "C" {
#include "espresso.h"
}

using namespace std;
using namespace tbb;
using namespace boost::multiprecision;
namespace mp = boost::multiprecision;
/*�򥻩w�q�P���c*/
#define ullg unsigned long long 

#define varnum 16 //�bbinary-input�����p�U�A�@��word�|��16��input
#define low cube.last_word[cube.output] // Last output word �̫�@��output��word
#define fow cube.first_word[cube.output] // First output word �Ĥ@��output��word
#define var_insert(set, e, value) (set[WHICH_WORD(e * 2)] |= (value+1) << WHICH_BIT(e * 2)) //�bpset����e���ܼƴ��Jvalue, 2 -> 11, 1 -> 01, 0 -> 10
#define var_remove(set, e) (set[WHICH_WORD(e * 2)] &= ~(3 << WHICH_BIT(e * 2))) //����pset����e���ܼ�

#define LastInputWord WHICH_WORD(NUMINPUTS * 2) //�̫�@�ӿ�J�Ҧb��word

#define Bat 500 //�Ω����ˬd���ˬd�I�A�Ω��ˬdSym-Table����ٱ��p sym�ϥ�


/*K-�S�x�Ȩϥ�*/
#define TABLE vector<ksig_value> //

struct ksig_value { //K-�S�x�Ȫ��V�q���
	vector<cpp_int> distx;
	ksig_value(int n) : distx(n) { }
	ksig_value() : distx(1) { } //default

	ksig_value& operator= (const ksig_value a) {
		if (distx.size() == a.distx.size()) {
			for (int i = 0; i < distx.size(); i++) {
				distx[i] = a.distx[i];
			}
		}
		return *this;
	}

	ksig_value operator+= (const ksig_value a) {
		if (distx.size() == a.distx.size()) {
			for (int i = 0; i < distx.size(); i++) {
				distx[i] += a.distx[i];
			}
		}
		return *this;
	}

	bool operator== (const ksig_value& a) const {
		if (distx.size() != a.distx.size())
			return false;

		for (int i = 0; i < distx.size(); i++) {
			if (distx[i] != a.distx[i]) {
				return false;
			}
		}
		return true;
	}
	bool operator!= (const ksig_value& a) const {
		if (distx.size() != a.distx.size())
			return true;

		for (int i = 0; i < distx.size(); i++) {
			if (distx[i] != a.distx[i]) {
				return true;
			}
		}
		return false;
	}
};

/*tool.cpp*/ char getvar(pset a, int position); //����pset��position���ȡA��X��?�B0�B1�B2(00�B10�B01�B11)
/*tool.cpp*/ string pstring(pset s, int n); //��Xpset�A�Ppbv1�A���L���׭���
/*tool.cpp*/ void sf_bm_printf(pset_family A); //��Xpset_family�A���ΤFpstring
/*tool.cpp*/ bool output_intersect(pset a, pset b); //�ˬd���pset��output�O�_�ۥ�
/*tool.cpp*/ pcover sep_sup_output(pset_family T, int i, pset sup); //�Ω�N��X���}
/*tool.cpp*/ void set_cube(int in, int out); //���s�]�wcube�����
/*tool.cpp*/ pset_family sf_rand(pset_family A); //�Npset_family�����ǥ���
/*tool.cpp*/ pset set_dc(pset r, pset a); //�^��pset��don't care(dc)��mask
/*tool.cpp*/ pset set_ndc(pset r, pset a); //�^��pset���Ddc��mask
/*tool.cpp*/ int dc_count(pset a); //�p��pset��dc���ƶq
/*tool.cpp*/ int dc_count(pset r, pset a); //�p��pset��dc���ƶq�A�æ^��dc��mask
/*tool.cpp*/ int set_ordc(pset a, pset b); //�p����pset���ܼƨ�@��dc���ƶq
/*tool.cpp*/ int set_anddc(pset a, pset b); //�p����pset���ܼƳ���dc���ƶq
/*tool.cpp*/ cpp_int Cxy(int x, int y); //�զX�ƾ�
/*tool.cpp*/ pset_family support_set(pset_family F, pset_family R); //�p���ƪ�support set
/*tool.cpp*/ void print_TABLE(TABLE t); //��XTABLE�榡
/*tool.cpp*/ pset_family sf_parallel_sort(pset_family pf); //�ϥ�TBB��@����Ƨ�pset_family�A��DC��J�ƶq���h��->���֪�

/*���P���󪺶Z���ˬd*/
/*tool.cpp*/ bool notcdist(pset a, pset b); //�p����pset�S���ۥ�^��1�A�Ϥ�0
/*tool.cpp*/ int cdist1(pset a, pset b); //�p����pset�Z����1�^��1�A�Ϥ�0
/*tool.cpp*/ int cdist123(pset a, pset b, int dist_on[2]); //�p����pset���Z���A�j��2���Z���Ҧ^��3�A�p�󵥩�2���Z���|�^�ǶZ�����ͪ��ܼƦ�m
/*tool.cpp*/ int cdist123(pset a, pset b); //�p����pset���Z���A�j��2���Z���^��3
/*tool.cpp*/ int cdistN(pset a, pset b); //�p����pset�����Z��

/*����ˬd*/
/*symmetry.cpp*/ pset set_rsp(pset r, pset p, pset q); //�p��ϦV���i������ٰt�諸pset
/*symmetry.cpp*/ bool redundant(pset a, pset b); //�ˬd��epset�O�_�|���������t��
/*symmetry.cpp*/ void find_Symmetry(pset_family A); //��X��ePLA����ٲզX=>�ƶq(��ٶ��X�ƶq)
/*symmetry.cpp*/ int EorNE(pset_family A, int eorne); //�ˬd��ePLA����E(eorne = 0)��NE(eorne = 1)���ƶq�Aeorne=2����X��̼ƶq

/*�������`��*/
/*symmetry.cpp*/ void sym_naive(pset_family F, pset_family R); //��t��k��@

/*C++�æ�귽�禡�w��@*/
/*cppsym.cpp*/ void cpp_sym(pset_family F, pset_family R, int t_); //��ٰt���u�ϥΤ@�Ӥ�����
/*cppsym.cpp*/ void cpp_thread(pset_family F, pset_family R, int t_); //��Ƥ��Φ��T�w����������ƶq(t_)
/*cppsym.cpp*/ void cpp_dcp(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_); //�ϥ�DCP���X��l�ާ@
/*cppsym.cpp*/ void cpp_dcp_ordering(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);
/*cppsym.cpp*/ void cpp_dcp_symrate(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);

/*testingcppsym.cpp*/ void sym_testing(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_); //�ϥ�DCP���X��l�ާ@�H�Ψϥ�weight�N��ƶ���J�u�@��C
/*testingcppsym.cpp*/ void sym_testing_ordering(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);
/*testingcppsym.cpp*/ void sym_testing_symrate(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);

/*Intel TBB��@*/
/*tbbsym.cpp*/ void tbb_sym(pset_family F, pset_family R, int cs_row, int cs_col);

/*/openMP��@*/
/*ompsym.cpp*/ void omp_sym(pset_family F, pset_family R);

/*ksig.cpp*/ void cal_ksignature(pPLA PLA, int Mode);
/*ksig.cpp*/ void ksio(pPLA PLA, int mode, int ioA);

//���ʹ���
/*GPLA.cpp*/ void general_setup(int input, int output); // ESPRESSO cube_setup
/*GPLA.cpp*/ void GeneralPLA_allrandom(pset_family& F, pset_family& R, int input_size, int onset_size, int offset_size, int on_off);
/*GPLA.cpp*/ bool GeneralPLA_nonSym(pset_family& F, pset_family& R, int input_size, int onset_size, int offset_size, int on_off);

/*�����פ夤���P������*/
/*experiment.cpp*/ void symmetry_benchmarking(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize); //����benchmarking���� (MCNC, ISCAS) 
/*experiment.cpp*/ void symmetry_LargeCircuit(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize); //���դj�����q������ (�ͦ����j���q��) (���O�ѼƤ������]�t -R)
/*experiment.cpp*/ void testing_chunksize(pset_family& F, pset_family& R, int ordering, int tsize); //�ϥΤ��Pchunksize��q���i�����
/*experiment.cpp*/ void testing_threadnum(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering); //�ϥΤ��P������ƶq��q���i�����
/*experiment.cpp*/ void testing_ordering(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int tsize); //�ϥΤ��P��ƶ���J�u�@��C���覡��q���i�����
/*experiment.cpp*/ void testing_symrate(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering); //�ϥγ�@������̧��ˬd��ƶ�����٩ʪ��U�����(�b�פ夤pdc�Pc3540_g409�q���ϥ�)