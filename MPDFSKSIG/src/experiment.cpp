#include "header.h"
#include "main.h"

//����benchmarking���� (MCNC, ISCAS) 
void symmetry_benchmarking(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering, int tsize) {
	pset_family F_ = sf_copy(sf_new(F->count, F->sf_size), F);
	pset_family R_ = sf_copy(sf_new(R->count, R->sf_size), R);

	sym_naive(F_, R_);
	tbb_sym(F, R, chunk_row, chunk_col);
	omp_sym(F, R);
	printf("unsort�q�� ");
	if (ordering == 2) {
		sym_testing(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	else {
		cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	F = sf_parallel_sort(F);
	R = sf_parallel_sort(R);
	printf("sort�q�� ");
	if (ordering == 2) {	
		sym_testing(F, R, chunk_row, chunk_col, ordering, tsize);
	}
	else {
		
		cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize);
	}

	FREE(F_);
	FREE(R_);
}

//���դj�����q������ (�ͦ����j���q��) (���O�ѼƤ������]�t -R)
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

//�ϥΤ��Pchunksize��q���i�����
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

//�ϥΤ��P������ƶq��q���i�����
void testing_threadnum(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering) {
	vector<int> tsize = { 1, 2, 4, 8, 16, 24 };

	for (int i = 0; i < tsize.size(); i++) {
		printf("������ƶq��%d ", tsize[i]);
		if (ordering == 2) {
			sym_testing(F, R, chunk_row, chunk_col, ordering, tsize[i]);
		}
		else {
			cpp_dcp(F, R, chunk_row, chunk_col, ordering, tsize[i]);
		}
	}
}

//�ϥΤ��P��ƶ���J�u�@��C���覡��q���i�����
void testing_ordering(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int tsize) {

	printf("�H��V�N��ƶ���J�u�@��C\n");
	cpp_dcp_ordering(F, R, chunk_row, chunk_col, 0, tsize);

	printf("�H�צV�N��ƶ���J�u�@��C\n");
	cpp_dcp_ordering(F, R, chunk_row, chunk_col, 1, tsize);

	printf("�H�v���N��ƶ���J�u�@��C\n");
	sym_testing_ordering(F, R, chunk_row, chunk_col, 2, tsize);
}

//�ϥγ�@������̧��ˬd��ƶ�����٩ʪ��U�����(�b�פ夤pdc�Pc3540_g409�q���ϥ�)
void testing_symrate(pset_family& F, pset_family& R, int chunk_row, int chunk_col, int ordering) {

	if (ordering == 2) {
		sym_testing_symrate(F, R, chunk_row, chunk_col, ordering, 1);
	}
	else {
		cpp_dcp_symrate(F, R, chunk_row, chunk_col, ordering, 1);
	}		
}
