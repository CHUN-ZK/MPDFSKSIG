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
/*基本定義與結構*/
#define ullg unsigned long long 

#define varnum 16 //在binary-input的情況下，一個word會有16個input
#define low cube.last_word[cube.output] // Last output word 最後一個output的word
#define fow cube.first_word[cube.output] // First output word 第一個output的word
#define var_insert(set, e, value) (set[WHICH_WORD(e * 2)] |= (value+1) << WHICH_BIT(e * 2)) //在pset中第e個變數插入value, 2 -> 11, 1 -> 01, 0 -> 10
#define var_remove(set, e) (set[WHICH_WORD(e * 2)] &= ~(3 << WHICH_BIT(e * 2))) //移除pset中第e個變數

#define LastInputWord WHICH_WORD(NUMINPUTS * 2) //最後一個輸入所在的word

#define Bat 500 //用於對稱檢查的檢查點，用於檢查Sym-Table的對稱情況 sym使用


/*K-特徵值使用*/
#define TABLE vector<ksig_value> //

struct ksig_value { //K-特徵值的向量表示
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

/*tool.cpp*/ char getvar(pset a, int position); //拿取pset中position的值，輸出為?、0、1、2(00、10、01、11)
/*tool.cpp*/ string pstring(pset s, int n); //輸出pset，同pbv1，但無長度限制
/*tool.cpp*/ void sf_bm_printf(pset_family A); //輸出pset_family，應用了pstring
/*tool.cpp*/ bool output_intersect(pset a, pset b); //檢查兩個pset的output是否相交
/*tool.cpp*/ pcover sep_sup_output(pset_family T, int i, pset sup); //用於將輸出分開
/*tool.cpp*/ void set_cube(int in, int out); //重新設定cube的資料
/*tool.cpp*/ pset_family sf_rand(pset_family A); //將pset_family的順序打亂
/*tool.cpp*/ pset set_dc(pset r, pset a); //回傳pset中don't care(dc)的mask
/*tool.cpp*/ pset set_ndc(pset r, pset a); //回傳pset中非dc的mask
/*tool.cpp*/ int dc_count(pset a); //計算pset中dc的數量
/*tool.cpp*/ int dc_count(pset r, pset a); //計算pset中dc的數量，並回傳dc的mask
/*tool.cpp*/ int set_ordc(pset a, pset b); //計算兩個pset中變數其一為dc的數量
/*tool.cpp*/ int set_anddc(pset a, pset b); //計算兩個pset中變數都為dc的數量
/*tool.cpp*/ cpp_int Cxy(int x, int y); //組合數學
/*tool.cpp*/ pset_family support_set(pset_family F, pset_family R); //計算函數的support set
/*tool.cpp*/ void print_TABLE(TABLE t); //輸出TABLE格式
/*tool.cpp*/ pset_family sf_parallel_sort(pset_family pf); //使用TBB實作平行排序pset_family，由DC輸入數量較多的->較少的

/*不同條件的距離檢查*/
/*tool.cpp*/ bool notcdist(pset a, pset b); //計算兩個pset沒有相交回傳1，反之0
/*tool.cpp*/ int cdist1(pset a, pset b); //計算兩個pset距離為1回傳1，反之0
/*tool.cpp*/ int cdist123(pset a, pset b, int dist_on[2]); //計算兩個pset的距離，大於2的距離皆回傳3，小於等於2的距離會回傳距離產生的變數位置
/*tool.cpp*/ int cdist123(pset a, pset b); //計算兩個pset的距離，大於2的距離回傳3
/*tool.cpp*/ int cdistN(pset a, pset b); //計算兩個pset間的距離

/*對稱檢查*/
/*symmetry.cpp*/ pset set_rsp(pset r, pset p, pset q); //計算反向的可移除對稱配對的pset
/*symmetry.cpp*/ bool redundant(pset a, pset b); //檢查當前pset是否會有移除的配對
/*symmetry.cpp*/ void find_Symmetry(pset_family A); //輸出當前PLA的對稱組合=>數量(對稱集合數量)
/*symmetry.cpp*/ int EorNE(pset_family A, int eorne); //檢查當前PLA中的E(eorne = 0)或NE(eorne = 1)的數量，eorne=2為輸出兩者數量

/*單執行緒循序*/
/*symmetry.cpp*/ void sym_naive(pset_family F, pset_family R); //原演算法實作

/*C++並行資源函式庫實作*/
/*cppsym.cpp*/ void cpp_sym(pset_family F, pset_family R, int t_); //對稱配對表只使用一個互斥鎖
/*cppsym.cpp*/ void cpp_thread(pset_family F, pset_family R, int t_); //資料切割成固定切成執行緒數量(t_)
/*cppsym.cpp*/ void cpp_dcp(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_); //使用DCP結合原子操作
/*cppsym.cpp*/ void cpp_dcp_ordering(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);
/*cppsym.cpp*/ void cpp_dcp_symrate(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);

/*testingcppsym.cpp*/ void sym_testing(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_); //使用DCP結合原子操作以及使用weight將資料塊放入工作佇列
/*testingcppsym.cpp*/ void sym_testing_ordering(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);
/*testingcppsym.cpp*/ void sym_testing_symrate(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_);

/*Intel TBB實作*/
/*tbbsym.cpp*/ void tbb_sym(pset_family F, pset_family R, int cs_row, int cs_col);

/*/openMP實作*/
/*ompsym.cpp*/ void omp_sym(pset_family F, pset_family R);

/*ksig.cpp*/ void cal_ksignature(pPLA PLA, int Mode);
/*ksig.cpp*/ void ksio(pPLA PLA, int mode, int ioA);

//產生測資
/*GPLA.cpp*/ void general_setup(int input, int output); // ESPRESSO cube_setup
/*GPLA.cpp*/ void GeneralPLA_allrandom(pset_family& F, pset_family& R, int input_size, int onset_size, int offset_size, int on_off);
/*GPLA.cpp*/ bool GeneralPLA_nonSym(pset_family& F, pset_family& R, int input_size, int onset_size, int offset_size, int on_off);

/*對應論文中不同的實驗*/
/*experiment.cpp*/ void symmetry_benchmarking(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize); //測試benchmarking實驗 (MCNC, ISCAS) 
/*experiment.cpp*/ void symmetry_LargeCircuit(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize); //測試大型的電路實驗 (生成的大型電路) (指令參數中必須包含 -R)
/*experiment.cpp*/ void testing_chunksize(pset_family& F, pset_family& R, int ordering, int tsize); //使用不同chunksize對電路進行實驗
/*experiment.cpp*/ void testing_threadnum(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering); //使用不同執行緒數量對電路進行實驗
/*experiment.cpp*/ void testing_ordering(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int tsize); //使用不同資料塊放入工作佇列的方式對電路進行實驗
/*experiment.cpp*/ void testing_symrate(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering); //使用單一執行緒依序檢查資料塊的對稱性的下降比例(在論文中pdc與c3540_g409電路使用)