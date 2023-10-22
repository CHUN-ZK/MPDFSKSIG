#include "oneapi/tbb/spin_mutex.h"
#include "header.h"

/*
	對稱性偵測(Intel TBB)
*/

class TBB_SYM {
private:
	pset_family F;			/*on-set積項*/
	pset_family R;			/*off-set積項*/
	pset_family sym_table;	/*對稱配對表*/
	pset S;					/*非對稱輸入集合*/
	int tsize;				/*執行緒數量*/
	int chunk_size;			/*資料塊size*/
	int chunk_row, chunk_col;
	vector<tbb::spin_rw_mutex> sym_lock;/*互斥鎖陣列*/
	bool is_not_sym;		/*用於表示對稱配對表是否為空集合*/
public:
	/*初始化*/
	TBB_SYM(pset_family F_, pset_family R_, int cs_row, int cs_col) : F(F_), R(R_), sym_lock(NUMINPUTS), is_not_sym(false), chunk_row(cs_row), chunk_col(cs_col) {
		S = set_full(NUMINPUTS * 2); 
		sym_table = sf_new(NUMINPUTS, NUMINPUTS * 2);
		for (int i = 0; i < NUMINPUTS; i++) {
			sf_addset(sym_table, S);
		}
		tsize = thread::hardware_concurrency();
		chunk_size = ceil((double)F->count / tsize);
	}

	void main_task() {
		clock_t start_time = clock();

		/*對稱檢查*/
		//static tbb::affinity_partitioner ap;
		tbb::parallel_for(tbb::blocked_range2d<int>(0, F->count, chunk_row, 0, R->count, chunk_col),
			[=](const tbb::blocked_range2d<int>& r) {
				pset k = set_full(NUMINPUTS * 2);
				int t = 0;
				bool b;
				for (int i = r.rows().begin(), i_end = r.rows().end(); i < i_end && !is_not_sym; i++) {
					if (!redundant(F->data + F->wsize * i, S)) {
						continue;
					}
					for (int j = r.cols().begin(), j_end = r.cols().end(); j < j_end; j++) {
						b = symmetry(i, j, k);


						if (b) {
							t++;
							if (t == Bat) {
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
				}
				FREE(k);
			}, tbb::auto_partitioner());

		allflip();
		printf("Intel TBB time = %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
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
					if (getvar(fp, dist_on[0]) == '1') {
						k = set_rsp(k, fp, rp);
					}
					else {
						k = set_rsp(k, rp, fp);
					}
					remove_pair1(k, dist_on[0]);
				}
				return true;
			}
		}
		return false;
	}

	/*lock all symmetry table and remove one line symmetry pair*/
	void remove_pair1(pset k, int var) {
		tbb::spin_rw_mutex::scoped_lock lock(sym_lock[var]);
		pset p = (sym_table->data + sym_table->wsize * var);
		set_and(p, p, k);
	}

	/*lock all symmetry table and remove E or NE pair*/
	void remove_pair2(int position, int var) {
		tbb::spin_rw_mutex::scoped_lock lock(sym_lock[var]);
		set_remove((sym_table->data + sym_table->wsize * var), position);
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
	~TBB_SYM() {
		sf_free(sym_table);
		FREE(S);
	}
};

void tbb_sym(pset_family F, pset_family R, int cs_row, int cs_col) {
	TBB_SYM sym(F, R, cs_row, cs_col);
	sym.main_task();
}
