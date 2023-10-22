#include <shared_mutex>
#include "header.h"

/*
	對稱性偵測(C++並行資源函式庫)
*/

class CPP_SYM {
protected:
	pset_family sym_table;		//對稱配對表
	pset_family F;				//on-set積項
	pset_family R;				//off-set積項
	pset S;						//非對稱輸入集合
	bool is_not_sym;			//用於表示對稱配對表是否為空集合
	int tsize;					//執行緒數量
private:
	mutable shared_mutex sym_lock;
	int chunk_size;				//依據執行緒數量計算資料塊大小

public:
	CPP_SYM(pset_family F_, pset_family R_, int t_) : F(F_), R(R_), is_not_sym(false), tsize(t_) {
		S = set_full(NUMINPUTS * 2); //store variable without symmetry
		sym_table = sf_new(NUMINPUTS, NUMINPUTS * 2);
		for (int i = 0; i < NUMINPUTS; i++) {
			sf_addset(sym_table, S);
		}
		chunk_size = ceil((double)F->count / tsize);
	}

	/*移除多組對稱配對*/
	void remove_pair1(pset k, int var) {
		lock_guard<shared_mutex> lock(sym_lock);
		pset p = (sym_table->data + sym_table->wsize * var);
		set_and(p, p, k);
	}

	/*移除一組對稱配對*/
	void remove_pair2(int position, int var) {
		lock_guard<shared_mutex> lock(sym_lock);
		set_remove((sym_table->data + sym_table->wsize * var), position);
	}

	void main_task() {
		clock_t start_time = clock();
		vector<thread> vec;
		
		/*對稱檢查*/
		for (int i = 0; i < tsize; i++) {
			vec.emplace_back(&CPP_SYM::symmetry, this, i);
		}
		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();
		printf("CPP + 1 Mutex time = %.3f\n", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}

	void symmetry(int cid) {
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int rnum, t = 0;

		for (int i = cid * chunk_size, j = 0; i < F->count && j < chunk_size && !is_not_sym; i++, j++) {
			fp = F->data + F->wsize * i;
			if (!redundant(fp, S)) {
				continue;
			}
			foreachi_set(R, rnum, rp) {
				if (output_intersect(fp, rp)) {
					int dist_on[2] = { -1, -1 };
					int dist = cdist123(fp, rp, dist_on);
					if (dist > 2) continue;
					else if (!is_in_set(S, dist_on[0] * 2)) {
						continue;
					}
					else if (dist == 2) {
						if (is_in_set(S, dist_on[1] * 2)) {
							int position = dist_on[1] * 2;
							if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
								position += 1;
							}
							remove_pair2(position, dist_on[0]);
						}
					}
					else if (dist == 1) {
						k = set_fill(k, NUMINPUTS * 2);
						if (getvar(fp, dist_on[0]) == '1') {
							k = set_rsp(k, fp, rp);
						}
						else {
							k = set_rsp(k, rp, fp);
						}
						remove_pair1(k, dist_on[0]);
					}
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
				for (int j = i + 1; j <
					NUMINPUTS; j++) {
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

	~CPP_SYM() {
		sf_free(sym_table);
		FREE(S);
	}
};

class CPP_THREAD : public CPP_SYM {
protected:
	vector<shared_mutex> sym_lock; //互斥鎖陣列
private:
	int chunk_size;

public:
	CPP_THREAD(pset_family F_, pset_family R_, int t_) : CPP_SYM(F_, R_, t_), sym_lock(NUMINPUTS) {
		chunk_size = ceil((double)F->count / tsize);
	}

	/*lock all symmetry table and remove one line symmetry pair*/
	void remove_pair1(pset k, int var) {
		lock_guard<shared_mutex> lock(sym_lock[var]);
		pset p = (sym_table->data + sym_table->wsize * var);
		set_and(p, p, k);
	}

	/*lock all symmetry table and remove E or NE pair*/
	void remove_pair2(int position, int var) {
		lock_guard<shared_mutex> lock(sym_lock[var]);
		set_remove((sym_table->data + sym_table->wsize * var), position);
	}

	void main_task() {
		clock_t start_time = clock();
		vector<thread> vec;
		for (int i = 0; i < tsize; i++) {
			vec.emplace_back(&CPP_THREAD::symmetry, this, i);
		}
		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();
		printf("CPP + Mutex Array time = %.3f\n", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}

	void symmetry(int cid) {
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int rnum, t = 0;
		for (int i = cid * chunk_size, j = 0; i < F->count && j < chunk_size && !is_not_sym; i++, j++) {
			fp = F->data + F->wsize * i;
			if (!redundant(fp, S)) {
				continue;
			}
			foreachi_set(R, rnum, rp) {
				if (output_intersect(fp, rp)) {
					int dist_on[2] = { -1, -1 };
					int dist = cdist123(fp, rp, dist_on);
					if (dist > 2) continue;
					else if (!is_in_set(S, dist_on[0] * 2)) {
						continue;
					}
					else if (dist == 2) {
						//is_in_set(S, dist_on[0] * 2)
						if (is_in_set(S, dist_on[1] * 2)) {
							int position = dist_on[1] * 2;
							if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
								position += 1;
							}
							remove_pair2(position, dist_on[0]);
						}
					}
					else if (dist == 1) {
						k = set_fill(k, NUMINPUTS * 2);
						if (getvar(fp, dist_on[0]) == '1') {
							k = set_rsp(k, fp, rp);
						}
						else {
							k = set_rsp(k, rp, fp);
						}
						remove_pair1(k, dist_on[0]);
					}
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
	}

	~CPP_THREAD() {
		sf_free(sym_table);
		FREE(S);
	}
};

class CPP_DCP : public CPP_THREAD {
protected:
	atomic<unsigned long long> chunk_index;
	int order_id, chunk_row, chunk_col;
	unsigned long long col, row;
	unsigned long long chunk_count;
	vector<pair<int, int>> chunk_table; //資料塊的Job Queue

public:
	CPP_DCP(pset_family F_, pset_family R_, int cs_row, int cs_col, int o, int t_) : CPP_THREAD(F_, R_, t_), chunk_row(cs_row), chunk_col(cs_col), order_id(o) {
		chunk_index.store(0);
		col = ceil((double)R->count / chunk_col);
		row = ceil((double)F->count / chunk_row);
		chunk_count = col * row;
		chunk_table.resize(chunk_count, make_pair<int, int>(-1, -1));
	}

	void main_task() {
		clock_t start_time = clock();
		
		/*資料塊放入Job Queue*/
		chunk_order(order_id);
		
		vector<thread> vec;
		int t;
		/*調整生成執行緒數量*/
		if (tsize * 2 > chunk_count) {
			t = ceil((double)chunk_count / 2);
		}
		else {
			t = tsize;
		}
		/*對稱檢查*/
		for (int i = 0; i < t; i++) {
			vec.emplace_back(&CPP_DCP::symmetry, this);
		}

		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();

		printf("Ours time =  %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}
	void symmetry() {
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int i, j, t1, t2, dist, times = 0;
		unsigned long long cid, chunk;
		int dist_on[2] = { -1, -1 };

		while (!is_not_sym && (cid = chunk_index.fetch_add(1)) < chunk_count) {
			i = chunk_table[cid].first * chunk_row;
			for (fp = F->data + F->wsize * i, t1 = 0; !is_not_sym && t1 < chunk_row && i < F->count; i++, t1++, fp += F->wsize) {
				if (!redundant(fp, S)) {
					continue;
				}
				j = chunk_table[cid].second * chunk_col;
				for (rp = R->data + R->wsize * j, t2 = 0; !is_not_sym && t2 < chunk_col && j < R->count; j++, t2++, rp += R->wsize) {
					if (output_intersect(fp, rp)) {
						dist = cdist123(fp, rp, dist_on);
						if (dist > 2) continue;
						else if (!is_in_set(S, dist_on[0] * 2)) {
							//skip2.fetch_add(1);
							continue;
						}
						else if (dist == 2) {
							if (is_in_set(S, dist_on[1] * 2)) {
								int position = dist_on[1] * 2;
								if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
									position += 1;
								}
								remove_pair2(position, dist_on[0]);
							}
						}
						else if (dist == 1) {
							k = set_fill(k, NUMINPUTS * 2);
							if (getvar(fp, dist_on[0]) == '1') {
								k = set_rsp(k, fp, rp);
							}
							else {
								k = set_rsp(k, rp, fp);
							}
							remove_pair1(k, dist_on[0]);
						}
						times++;
						if (times  == Bat) {
							times = 0;
							flip();
							
							if (sym_check()) {
								is_not_sym = true;
								return;
							}
						}
					}
				}
			}
		}
		FREE(k);
	}

	void main_ordering() {
		clock_t start_time = clock();

		/*資料塊放入Job Queue*/
		chunk_order(order_id);

		/*輸出資料塊放入工作佇列的時間*/
		printf("Ordering time =  %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);

		start_time = clock();

		vector<thread> vec;
		int t;
		/*調整生成執行緒數量*/
		if (tsize * 2 > chunk_count) {
			t = ceil((double)chunk_count / 2);
		}
		else {
			t = tsize;
		}
		/*對稱檢查*/
		for (int i = 0; i < t; i++) {
			vec.emplace_back(&CPP_DCP::symmetry, this);
		}

		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();

		printf("Ours time =  %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}

	void symrate() {
		chunk_order(order_id);
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int i, j, t1, t2, dist, times = 0;
		unsigned long long cid, chunk;
		int dist_on[2] = { -1, -1 };

		int total_sym = EorNE(sym_table, 2);
		printf("%.3f % \n", (double)total_sym / total_sym * 100);

		while (!is_not_sym && (cid = chunk_index.fetch_add(1)) < chunk_count) {
			i = chunk_table[cid].first * chunk_row;
			for (fp = F->data + F->wsize * i, t1 = 0; !is_not_sym && t1 < chunk_row && i < F->count; i++, t1++, fp += F->wsize) {
				if (!redundant(fp, S)) {
					continue;
				}
				j = chunk_table[cid].second * chunk_col;
				for (rp = R->data + R->wsize * j, t2 = 0; !is_not_sym && t2 < chunk_col && j < R->count; j++, t2++, rp += R->wsize) {
					if (output_intersect(fp, rp)) {
						dist = cdist123(fp, rp, dist_on);
						if (dist > 2) continue;
						else if (!is_in_set(S, dist_on[0] * 2)) {
							//skip2.fetch_add(1);
							continue;
						}
						else if (dist == 2) {
							if (is_in_set(S, dist_on[1] * 2)) {
								int position = dist_on[1] * 2;
								if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
									position += 1;
								}
								remove_pair2(position, dist_on[0]);
							}
						}
						else if (dist == 1) {
							k = set_fill(k, NUMINPUTS * 2);
							if (getvar(fp, dist_on[0]) == '1') {
								k = set_rsp(k, fp, rp);
							}
							else {
								k = set_rsp(k, rp, fp);
							}
							remove_pair1(k, dist_on[0]);
						}
						times++;
						if (times == Bat) {
							times = 0;
							flip();

							if (sym_check()) {
								is_not_sym = true;
								allflip();
								printf("%.3f % \n", (double)EorNE(sym_table, 2) / total_sym * 100);
								return;
							}
						}
					}
				}
			}
			allflip();
			printf("%.3f % \n", (double)EorNE(sym_table, 2) / total_sym * 100);
		}
		allflip();
		printf("%.3f % \n", (double)EorNE(sym_table, 2) / total_sym * 100);

		FREE(k);
	}

	/*chunksize = cs_row*cs_col、工作佇列內資料塊順序(橫向(0)、斜向(1))、最大執行緒數量為t_*/
	void chunk_order(int order) {
		if (order == 1) {
			int tmp = 0;
			for (int t = 0; t < col; t++) {
				for (int i = 0, j = (t - i); i < row && j >= 0 && (i + j) == t; i++, j--) {
					chunk_table[tmp].first = i;
					chunk_table[tmp].second = j;
					tmp++;
				}
			}
			for (int t = col; t <= row + col - 2; t++) {
				for (int i = t - (col - 1), j = (t - i); i < row && j >= 0 && (i + j) == t; i++, j--) {
					chunk_table[tmp].first = i;
					chunk_table[tmp].second = j;
					tmp++;
				}
			}
		}
		else {
			int tmp = 0;
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {
					chunk_table[tmp].first = i;
					chunk_table[tmp].second = j;
					tmp++;
				}
			}
		}
	}	
	~CPP_DCP() {
		sf_free(sym_table);
		FREE(S);
	}
};

/*只有一個互斥鎖*/
void cpp_sym(pset_family F, pset_family R, int t_) {
	CPP_SYM sym(F, R, t_);
	sym.main_task();
}

/*使用全部可用執行緒數量*/
void cpp_thread(pset_family F, pset_family R, int t_) {
	CPP_THREAD sym(F, R, t_);
	sym.main_task();
}

/*DCP*/
void cpp_dcp(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_) {
	CPP_DCP sym(F, R, cs_row, cs_col, o, t_);
	sym.main_task();
}

void cpp_dcp_ordering(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_) {
	CPP_DCP sym(F, R, cs_row, cs_col, o, t_);
	sym.main_ordering();
}

void cpp_dcp_symrate(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_) {
	CPP_DCP sym(F, R, cs_row, cs_col, o, t_);
	sym.symrate();
}