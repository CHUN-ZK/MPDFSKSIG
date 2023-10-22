#include "header.h"

class OMP_SYM {
private:
	pset_family F;			/*on-set積項*/
	pset_family R;			/*off-set積項*/
	pset_family sym_table;	/*對稱配對表*/
	pset S;					/*非對稱輸入集合*/
	int tsize;				/*執行緒數量*/
	vector < omp_lock_t> sym_lock;/*互斥鎖陣列*/
	bool is_not_sym;		/*用於表示對稱配對表是否為空集合*/

public:
	OMP_SYM(pset_family F_, pset_family R_) : F(F_), R(R_), sym_lock(NUMINPUTS), is_not_sym(false) {
		S = set_full(NUMINPUTS * 2); //store variable without symmetry
		sym_table = sf_new(NUMINPUTS, NUMINPUTS * 2);
		for (int i = 0; i < NUMINPUTS; i++) {
			sf_addset(sym_table, S);
		}
		tsize = thread::hardware_concurrency();
		for (int i = 0; i < NUMINPUTS; i++) {
			omp_init_lock(&sym_lock[i]);
		}
	}

	void main_task() {
		clock_t start_time = clock();
		int i, j, t;

		/*對稱檢查*/
		omp_set_nested(1);
#pragma omp parallel for  private(i,j,t) schedule(dynamic)
		for (i = 0; i < F->count; i++)
		{
			pset k = set_full(NUMINPUTS * 2);
			bool b;
			t = 0;
			for (j = 0; j < R->count; j++)
			{
				if (!redundant(F->data + F->wsize * i, S)) {
					continue;
				}
				b = symmetry(i, j, k);
				if (b) {
					t++;
					if (t >= Bat) {
						t = 0;
						flip();
						/*Stop if there is no symmetry*/
						if (sym_check()) {
							is_not_sym = true;			
							break;
						}
					}
				}
			}
			FREE(k);
		}
		allflip();
		printf("OpenMP time = %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}
	bool symmetry(int i, int j, pset k) {
		pset fp = F->data + F->wsize * i;
		pset rp = R->data + R->wsize * j;
		if (output_intersect(fp, rp)) {
			int dist_on[2] = { -1, -1 };
			int dist = cdist123(fp, rp, dist_on);
			if (is_in_set(S, dist_on[0] * 2) && dist < 3) {
				if (dist == 2) {
					if (is_in_set(S, dist_on[1] * 2)) {
						int position = dist_on[1] * 2;
						if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
							position += 1;
						}
						remove_pair2(position, dist_on[0]);
					}
				}
				else {
					//pset k = set_full(NUMINPUTS * 2);
					if (getvar(fp, dist_on[0]) == '1') {
						k = set_rsp(k, fp, rp);
					}
					else {
						k = set_rsp(k, rp, fp);
					}
					remove_pair1(k, dist_on[0]);
					//FREE(k);
				}
				return true;
			}
		}
		return false;
	}

	/*lock all symmetry table and remove one line symmetry pair*/
	void remove_pair1(pset k, int var) {
		omp_set_lock(&sym_lock[var]);
		pset p = (sym_table->data + sym_table->wsize * var);
		set_and(p, p, k);
		omp_unset_lock(&sym_lock[var]);
	}

	/*lock all symmetry table and remove E or NE pair*/
	void remove_pair2(int position, int var) {
		omp_set_lock(&sym_lock[var]);
		set_remove((sym_table->data + sym_table->wsize * var), position);
		omp_unset_lock(&sym_lock[var]);
	}

	/*將對稱配對表中的下矩陣統整到上矩陣*/
	void flip() {
		pset p, q;
		int i;

		foreachi_set(sym_table, i, p) {
			q = sym_table->data + i * sym_table->wsize;
			for (int j = i + 1; j < NUMINPUTS; j++) {
				q += sym_table->wsize;
				if (!is_in_set(q, i * 2) && is_in_set(p, j * 2)) {
					remove_pair2(j * 2, i);
				}
				if (!is_in_set(q, i * 2 + 1) && is_in_set(p, j * 2 + 1)) {
					remove_pair2(j * 2 + 1, i);
				}
			}
		}
	}

	/*將對稱配對表中上矩陣與下矩陣整合*/
	void allflip() {
		pset p, q;
		int i;

		foreachi_set(sym_table, i, p) {
			q = sym_table->data + i * sym_table->wsize;
			for (int j = i + 1; j < NUMINPUTS; j++) {
				q += sym_table->wsize;
				if (!is_in_set(q, i * 2)) {
					set_remove(p, j * 2);
				}
				if (!is_in_set(q, i * 2 + 1)) {
					set_remove(p, j * 2 + 1);
				}
				if (!is_in_set(p, j * 2)) {
					set_remove(q, i * 2);
				}
				if (!is_in_set(p, j * 2 + 1)) {
					set_remove(q, i * 2 + 1);
				}
			}
		}
	}

	/*檢查對稱配對表是否為空集合*/
	bool sym_check() {
		bool terminal = true;
		vector<bool> check(NUMINPUTS, false);
		pset p;

		for (int i = 0; i < NUMINPUTS; i++) {
			if (is_in_set(S, i * 2)) {
				p = sym_table->data + sym_table->wsize * i;
				for (int j = i + 1; j < NUMINPUTS; j++) {
					if (getvar(p, j) != '?') {
						check[i] = true;
						check[j] = true;
						terminal = false;
					}
				}
				if (!check[i]) {
					var_remove(S, i);
				}
			}
		}
		return terminal;
	}

	~OMP_SYM() {
		sf_free(sym_table);
		FREE(S);
		for (int i = 0; i < NUMINPUTS; i++) {
			omp_destroy_lock(&sym_lock[i]);
		}
	}
};

void omp_sym(pset_family F, pset_family R) {
	OMP_SYM sym(F, R);
	sym.main_task();
}
/*OpenMP版本計算KSO*/
/*C++並行函式庫版本*/
/*C++並行函式庫版本+DCP計算KSO與KSI*/