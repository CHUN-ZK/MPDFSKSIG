#include "header.h"
#include "main.h"

//測試benchmarking實驗 (MCNC, ISCAS) 
void symmetry_benchmarking(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize) {
	pset_family F_ = sf_copy(sf_new(F->count, F->sf_size), F);
	pset_family R_ = sf_copy(sf_new(R->count, R->sf_size), R);

	sym_naive(F_, R_);
	tbb_sym(F, R, chunk_row, chunk_col);
	omp_sym(F, R);
	printf("unsort電路 ");
	if (ordering == 2) {
		sym_testing(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	else {
		cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	F = sf_parallel_sort(F);
	R = sf_parallel_sort(R);
	printf("sort電路 ");
	if (ordering == 2) {	
		sym_testing(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	else {
		
		cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize);
	}

	FREE(F_);
	FREE(R_);
}

//測試大型的電路實驗 (生成的大型電路) (指令參數中必須包含 -R)
void symmetry_LargeCircuit(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize) {
	pset_family F_ = sf_copy(sf_new(F->count, F->sf_size), F);
	pset_family R_ = sf_copy(sf_new(R->count, R->sf_size), R);

	sym_naive(F_, R_);
	tbb_sym(F, R, chunk_row, chunk_col);
	omp_sym(F, R);
	if (ordering == 2) {
		sym_testing(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	else {
		cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize);
	}
}

//使用不同chunksize對電路進行實驗
void testing_chunksize(pset_family& F, pset_family& R, int ordering, int tsize) {
	vector<int> chunk_row = { 50, 200, 1000, 5000, 25000, 50, 200, 1000, 5000, 25000, (int)ceil((double)F->count / 10) };
	vector<int> chunk_col = { 50, 200, 1000, 5000, 25000, R->count, R->count, R->count,R->count, R->count, (int)ceil((double)R->count / 10) };

	for (int i = 0; i < chunk_row.size(); i++) {
		printf("chunksize = %d X %d ", chunk_row[i], chunk_col[i]);
		if (ordering == 2) {
			sym_testing(F, R, chunk_row[i], chunk_col[i], ordering, tsize);
		}
		else {
			cpp_dcp(F, R, chunk_row[i], chunk_col[i], ordering, tsize);
		}
	}
}

//使用不同執行緒數量對電路進行實驗
void testing_threadnum(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering) {
	vector<int> tsize = { 1, 2, 4, 8, 16, 24 };

	for (int i = 0; i < tsize.size(); i++) {
		printf("執行緒數量為%d ", tsize[i]);
		if (ordering == 2) {
			sym_testing(F, R, chunk_row, chunk_col, ordering, tsize[i]);
		}
		else {
			cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize[i]);
		}
	}
}

//使用不同資料塊放入工作佇列的方式對電路進行實驗
void testing_ordering(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int tsize) {

	printf("以橫向將資料塊放入工作佇列\n");
	cpp_dcp_ordering(F, R, chunk_row, chunk_col, 0, tsize);

	printf("以斜向將資料塊放入工作佇列\n");
	cpp_dcp_ordering(F, R, chunk_row, chunk_col, 1, tsize);

	printf("以權重將資料塊放入工作佇列\n");
	sym_testing_ordering(F, R, chunk_row, chunk_col, 2, tsize);
}

//使用單一執行緒依序檢查資料塊的對稱性的下降比例(在論文中pdc與c3540_g409電路使用)
void testing_symrate(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering) {

	if (ordering == 2) {
		sym_testing_symrate(F, R, chunk_row, chunk_col, ordering, 1);
	}
	else {
		cpp_dcp_symrate(F, R, chunk_row, chunk_col, ordering, 1);
	}		
}
