#include <shared_mutex>
#include <algorithm>
#include "header.h"

/*
	對稱性偵測(C++並行資源函式庫)
	使用權重(w1=資料塊內移除對稱配對數量, w2=資料塊內需進行移除對稱配對的輸入集合)
*/

struct chunk {
	pair<int, int> pid; //用於儲存區塊的座標
	int w1; //第一權重
	int w2; //第二權重
};

class SYM_TEST {
protected:
	pset_family sym_table;
	pset_family F;
	pset_family R;
	pset S;
	bool is_not_sym;
	int tsize, chunk_size;
	mutex mtx;
	atomic<unsigned long long> chunk_index;
	int order_id, chunk_row, chunk_col;
	unsigned long long col, row;
	unsigned long long chunk_count;
	atomic<int> runchunk;

private:
	vector<shared_mutex> sym_lock;
	vector<pair<int,int>> chunk_table;
	vector<chunk> chunk_table2;

public:
	SYM_TEST(pset_family F_, pset_family R_, int cs_row, int cs_col, int o, int t_) : F(F_), R(R_), is_not_sym(false), tsize(t_), sym_lock(NUMINPUTS), chunk_row(cs_row), chunk_col(cs_col), order_id(o) {
		S = set_full(NUMINPUTS * 2); //store variable without symmetry
		sym_table = sf_new(NUMINPUTS, NUMINPUTS * 2);
		for (int i = 0; i < NUMINPUTS; i++) {
			sf_addset(sym_table, S);
		}
		
		chunk_index.store(0);
		runchunk.store(0);
		col = ceil((double)R->count / chunk_col);
		row = ceil((double)F->count / chunk_row);
		chunk_count = col * row;
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

	void main_task2() {
		/*資料塊放入Job Queue*/
		chunk_table2.resize(chunk_count);
		for (int i = 0; i < chunk_count; i++) {
			chunk_table2[i].pid = make_pair<int, int>(-1, -1);
		}
		chunk_order(2);

		clock_t start_time = clock();

		chunk_index = 0;

		/*設定執行緒數量*/
		int t;
		if (tsize * 2 > chunk_count) {
			t = ceil((double)chunk_count / 2);
		}
		else {
			t = tsize;
		}
		
		/*對稱檢查*/
		vector<thread> vec;
		for (int i = 0; i < t; i++) {
			vec.emplace_back(&SYM_TEST::symmetry2, this);
		}

		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();
		printf("Ours(Weight sort) = %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table);
	}
	void main_task4() {
		/*資料塊放入Job Queue*/
		chunk_table2.resize(chunk_count);
		for (int i = 0; i < chunk_count; i++) {
			chunk_table2[i].pid = make_pair<int, int>(-1, -1);
		}
		chunk_order(4);

		chunk_index = 0;

		clock_t start_time = clock();

		/*設定執行緒數量*/
		int t;
		if (tsize * 2 > chunk_count) {
			t = ceil((double)chunk_count / 2);
		}
		else {
			t = tsize;
		}

		/*對稱檢查*/
		vector<thread> vec;
		for (int i = 0; i < t; i++) {
			vec.emplace_back(&SYM_TEST::symmetry4, this);
		}

		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();
		printf("Ours(前處理) = %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}

	void symmetry2() {
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int i, j, t1, t2, dist, times = 0;
		unsigned long long cid, chunk;
		int dist_on[2] = { -1, -1 };
		int run = 0;
		while (!is_not_sym && (cid = chunk_index.fetch_add(1)) < chunk_count) {
			i = chunk_table2[cid].pid.first * chunk_row;
			for (fp = F->data + F->wsize * i, t1 = 0; !is_not_sym && t1 < chunk_row && i < F->count; i++, t1++, fp += F->wsize) {
				if (!redundant(fp, S)) {
					//skip1.fetch_add(1);
					continue;
				}
				j = chunk_table2[cid].pid.second * chunk_col;
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
							//t.store(0);
							times = 0;
							flip();
							/*Stop if there is no symmetry*/
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
	void symmetry4() {
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int i, j, t1, t2, dist, times = 0;
		unsigned long long cid, chunk;
		int dist_on[2] = { -1, -1 };
		int run = 0;
		while (!is_not_sym && (cid = chunk_index.fetch_add(1)) < chunk_count) {
			if (chunk_table2[cid].w1 == 0)
				continue;
			i = chunk_table2[cid].pid.first * chunk_row;
			for (fp = F->data + F->wsize * i, t1 = 0; !is_not_sym && t1 < chunk_row && i < F->count; i++, t1++, fp += F->wsize) {
				if (!redundant(fp, S)) {
					continue;
				}
				j = chunk_table2[cid].pid.second * chunk_col;
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
							//t.store(0);
							times = 0;
							flip();
							/*Stop if there is no symmetry*/
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
		chunk_table2.resize(chunk_count);
		for (int i = 0; i < chunk_count; i++) {
			chunk_table2[i].pid = make_pair<int, int>(-1, -1);
		}
		chunk_order(2);

		chunk_index = 0;
		/*輸出資料塊放入工作佇列的時間*/
		printf("Ordering time =  %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);

		start_time = clock();

		/*設定執行緒數量*/
		int t;
		if (tsize * 2 > chunk_count) {
			t = ceil((double)chunk_count / 2);
		}
		else {
			t = tsize;
		}

		/*對稱檢查*/
		vector<thread> vec;
		for (int i = 0; i < t; i++) {
			vec.emplace_back(&SYM_TEST::symmetry2, this);
		}

		for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		allflip();
		printf("Ours(Weight sort) = %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table);
	}

	void symrate() {
		/*資料塊放入Job Queue*/
		chunk_table2.resize(chunk_count);
		for (int i = 0; i < chunk_count; i++) {
			chunk_table2[i].pid = make_pair<int, int>(-1, -1);
		}
		chunk_order(2);

		chunk_index = 0;
			
		pset fp, rp, k = set_new(NUMINPUTS * 2);
		int i, j, t1, t2, dist, times = 0;
		unsigned long long cid, chunk;
		int dist_on[2] = { -1, -1 };
		int run = 0;

		int total_sym = EorNE(sym_table, 2);
		printf("%.3f % \n", (double)total_sym / total_sym * 100);

		while (!is_not_sym && (cid = chunk_index.fetch_add(1)) < chunk_count) {
			i = chunk_table2[cid].pid.first * chunk_row;
			for (fp = F->data + F->wsize * i, t1 = 0; !is_not_sym && t1 < chunk_row && i < F->count; i++, t1++, fp += F->wsize) {
				if (!redundant(fp, S)) {
					continue;
				}
				j = chunk_table2[cid].pid.second * chunk_col;
				for (rp = R->data + R->wsize * j, t2 = 0; !is_not_sym && t2 < chunk_col && j < R->count; j++, t2++, rp += R->wsize) {
					if (output_intersect(fp, rp)) {
						dist = cdist123(fp, rp, dist_on);
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
						times++;
						if (times == Bat) {
							//t.store(0);
							times = 0;
							flip();
							/*Stop if there is no symmetry*/
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

	void chunk_weight() {
		pset fp, rp;
		int i, j, t1, t2, dist, times = 0;
		unsigned long long cid, chunk;
		int dist_on[2] = { -1, -1 };
		int d1, d2, w2 = 0;
		pset weight2 = set_new(NUMINPUTS);
		while ((cid = chunk_index.fetch_add(1)) < chunk_count) {
			d1 = 0;
			d2 = 0;
			w2 = 0;
			set_clear(weight2, NUMINPUTS);

			i = chunk_table2[cid].pid.first * chunk_row;
			for (fp = F->data + F->wsize * i, t1 = 0; !is_not_sym && t1 < chunk_row && i < F->count; i++, t1++, fp += F->wsize) {
				j = chunk_table2[cid].pid.second * chunk_col;
				for (rp = R->data + R->wsize * j, t2 = 0; !is_not_sym && t2 < chunk_col && j < R->count; j++, t2++, rp += R->wsize) {
					if (output_intersect(fp, rp)) {
						dist = cdist123(fp, rp, dist_on);
						if (dist == 2) {
							d2++;
							set_insert(weight2, dist_on[0]);
						}
						else if (dist == 1) {
							d1++;
							set_insert(weight2, dist_on[0]);
						}
					}
				}
			}
			for (i = 0; i < NUMINPUTS; i++) {
				if (is_in_set(weight2, i)) {
					w2++;
				}
			}
			if (d1 * 2 + d2 > 0) {
				runchunk.fetch_add(1);
				chunk_table2[cid].w1 = d1 * 2 + d2;
				chunk_table2[cid].w2 = w2;
			}
		}
		FREE(weight2);
	}

	void chunk_order(int order) {
		if (order == 2) {
			int tmp = 0;
			for (int r = 0; r < row; r++) {
				for (int c = 0; c < col; c++) {
					chunk_table2[tmp].pid.first = r;
					chunk_table2[tmp].pid.second = c;
					tmp++;
				}
			}

			int t;
			if (tsize * 2 > chunk_count) {
				t = ceil((double)chunk_count / 2);
			}
			else {
				t = tsize;
			}
			vector<thread> vec;
			for (int i = 0; i < t; i++) {
				vec.emplace_back(&SYM_TEST::chunk_weight, this);
			}

			for_each(vec.begin(), vec.end(), mem_fn(&thread::join));

			chunk_count = runchunk;
			sort(chunk_table2.begin(), chunk_table2.end(), [](chunk a, chunk b) {
				if (a.w1 == b.w1) {
					return a.w2 > b.w2;
				}
				else {
					return a.w1 > b.w1;
				}
				});
		}
		else if (order == 4) {
			int tmp = 0;
			for (int t = 0; t < col; t++) {
				for (int i = 0, j = (t - i); i < row && j >= 0 && (i + j) == t; i++, j--) {
					chunk_table2[tmp].pid.first = i;
					chunk_table2[tmp].pid.second = j;
					tmp++;
				}
			}
			for (int t = col; t <= row + col - 2; t++) {
				for (int i = t - (col - 1), j = (t - i); i < row && j >= 0 && (i + j) == t; i++, j--) {
					chunk_table2[tmp].pid.first = i;
					chunk_table2[tmp].pid.second = j;
					tmp++;
				}
			}
			int t;
			if (tsize * 2 > chunk_count) {
				t = ceil((double)chunk_count / 2);
			}
			else {
				t = tsize;
			}
			vector<thread> vec;
			for (int i = 0; i < t; i++) {
				vec.emplace_back(&SYM_TEST::chunk_weight, this);
			}

			for_each(vec.begin(), vec.end(), mem_fn(&thread::join));
		}
		else {
			 std::cout << "nothing" << std::endl;
		}
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
					//set_remove(p, j * 2);
					remove_pair2(j * 2, i);
				}
				if (!is_in_set(q, i * 2 + 1) && is_in_set(p, j * 2 + 1)) {
					//set_remove(p, j * 2 + 1);
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

	~SYM_TEST() {
		sf_free(sym_table);
		FREE(S);
	}
};

void sym_testing(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_) {
	SYM_TEST sym(F, R, cs_row, cs_col, o, t_);
	if (o == 2)
		sym.main_task2();
	else if (o == 4)
		sym.main_task4();
}

void sym_testing_ordering(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_) {
	SYM_TEST sym(F, R, cs_row, cs_col, o, t_);
	if (o == 2)
		sym.main_ordering();
}

void sym_testing_symrate(pset_family F, pset_family R, int cs_row, int cs_col, int o, int t_) {
	SYM_TEST sym(F, R, cs_row, cs_col, o, t_);
	if (o == 2)
		sym.symrate();
}