#include "header.h"
#include "main.h"
/*
	對稱性偵測使用函式
	對稱性偵測(單執行緒循序)
*/

/*計算可移除的對稱配對*/
pset set_rsp(pset r, pset p, pset q) {
	unsigned int p_dc, q_dc, ABC, C, AB;
	int tmp;

	if ((tmp = cube.inword) != -1) {
		p_dc = (p[tmp] >> 1) & p[tmp] & cube.inmask;
		p_dc = p_dc | (p_dc << 1);
		q_dc = (q[tmp] >> 1) & q[tmp] & cube.inmask;
		q_dc = q_dc | (q_dc << 1);
		ABC = p_dc | q_dc;
		C = ABC & (~q_dc);
		AB = ABC & (~C);

		r[tmp] = ~((AB & p[tmp]) | (~q[tmp] & C));
	}

	for (int i = 1; i < tmp; i++) {
		p_dc = (p[i] >> 1) & p[i] & DISJOINT;
		p_dc = p_dc | (p_dc << 1);
		q_dc = (q[i] >> 1) & q[i] & DISJOINT;
		q_dc = q_dc | (q_dc << 1);
		ABC = p_dc | q_dc;
		C = ABC & (~q_dc);
		AB = ABC & (~C);

		r[i] = ~((AB & p[i]) | (~q[i] & C));
	}

	return r;
}

/*檢查積項a是否為多餘積項*/
bool redundant(pset a, pset S) {
	unsigned int x;
	int tmp;

	if ((tmp = cube.inword) != -1) {
		x = (a[tmp] >> 1) ^ a[tmp]; //取非dc部分
		x = x & S[tmp] & cube.inmask;
		if (x)
			return true;
	}
	for (int i = 1; i < tmp; i++) {
		x = (a[i] >> 1) ^ a[i];
		x = x & S[i] & DISJOINT;
		if (x)
			return true;
	}
	return false;
}

/*輸出對稱表A中對稱的組合=>數量(對稱數量)*/
void find_Symmetry(pset_family Sym) {
	vector<int> times(NUMINPUTS + 1, 0);
	vector<int> isSym(NUMINPUTS, 1);
	int cnt = 0, count;
	for (int i = 0; i < NUMINPUTS; i++) {
		count = 0;
		pset p = Sym->data + Sym->wsize * i;
		for (int j = i; j < NUMINPUTS; j++) {
			if (isSym[j]) {
				if (is_in_set(p, j * 2)) {
					isSym[j] = 0;
					count++;
				}
				else if (is_in_set(p, j * 2 + 1)) {
					isSym[j] = 0;
					count++;
				}
			}
		}
		if (count != 0)
			times[count]++;
	}

	if (verbose) {
		printf("Sym = ");
		for (int i = 0; i <= NUMINPUTS; i++) {
			if (times[i] != 0) {
				printf("%d(%d)", times[i], i);
			}
		}
		printf("\n");
	}
}

/*檢查對稱配對表A中的對稱的數量，0 -> NE , 1 -> E, 2 -> E & NE*/
int EorNE(pset_family A, int eorne) {
	int i, j, cnt = 0;
	pset p;
	if (eorne == 2) {
		foreachi_set(A, i, p) {
			for (j = i + 1; j < NUMINPUTS; j++) {
				if (is_in_set(p, j * 2)) {
					cnt++;
				}
				if (is_in_set(p, j * 2 + 1)) {
					cnt++;
				}
			}
		}
	}
	else {
		foreachi_set(A, i, p) {
			for (j = i + 1; j < NUMINPUTS; j++) {
				if (is_in_set(p, j * 2 + eorne)) {
					cnt++;
				}
			}
		}
	}
	return cnt;
}

/*單執行緒循序方法*/
class Sym {
protected:
	pset_family F;  /*on-set*/
	pset_family R;	/*off-set*/
	pset_family sym_table; /*對稱配對表*/
	pset S; /*非對稱輸入集合*/
	bool state_is_change; /*用於檢查對稱配對表是否有更改*/

public:
	/*初始化*/
	Sym(pset_family F_, pset_family R_) : F(F_), R(R_) {
		S = set_full(NUMINPUTS * 2);
		sym_table = sf_new(NUMINPUTS, NUMINPUTS * 2);
		for (int i = 0; i < NUMINPUTS; i++) {
			sf_addset(sym_table, S);
		}
		state_is_change = 0;
	}

	void main_naive() {
		clock_t start_time = clock();
		/*對稱檢查*/
		symmetry_naive();
		printf("Naive time = %.3f ", (clock() - start_time) / (double)CLOCKS_PER_SEC);
		find_Symmetry(sym_table); //n(m) -> 有n個對稱群組的成員數量為m
	}

	void symmetry_naive() {
		pset fp, rp, p, k = set_full(NUMINPUTS * 2);;
		int t = 0, fnum, rnum;
		foreachi_set(F, fnum, fp) {
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
						int position = dist_on[1] * 2, position2 = dist_on[0] * 2;
						if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
							position += 1;
							position2 += 1;
						}
						p = sym_table->data + sym_table->wsize * dist_on[0];
						set_remove(p, position);

						p = sym_table->data + sym_table->wsize * dist_on[1];
						set_remove(p, position2);
					}
					else if (dist == 1) {
						k = set_fill(k, NUMINPUTS * 2);
						if (getvar(fp, dist_on[0]) == '1') {
							k = set_rsp(k, fp, rp);
						}
						else {
							k = set_rsp(k, rp, fp);
						}
						int i;
						for (i = 0, p = sym_table->data; i < NUMINPUTS; i++, p += sym_table->wsize) {
							if (i != dist_on[0]) {
								if (!is_in_set(k, i * 2)) {
									set_remove(p, dist_on[0] * 2);
								}
								if (!is_in_set(k, i * 2 + 1)) {
									set_remove(p, dist_on[0] * 2 + 1);
								}
							}
						}
						p = sym_table->data + sym_table->wsize * dist_on[0];
						set_and(p, p, k);
					}
					
					if (sym_check()) {
						return;
					}
				}
			}

			if (sym_check_naive()) {
				return;
			}
			/*Reduction*/
			if (state_is_change) {
				foreachi_set(R, rnum, rp) {
					if (!redundant(rp, S)) {
						rp -= R->wsize;
						sf_delset(R, rnum);
						rnum--;
					}
				}
			}
		}
		FREE(k);
	}

	/*檢查對稱配對表是否為空集合*/
	bool sym_check() {
		bool terminal = true, nosym;
		vector<bool> check(NUMINPUTS, false);
		pset p;
		int i;

		for (i = 0, p = sym_table->data; terminal && i < NUMINPUTS; i++, p += sym_table->wsize) {
			nosym = setp_equal(p, cube.var_mask[i]);

			if (!nosym) {
				terminal = false;
			}
		}
		return terminal;
	}

	/*檢查對稱配對表是否為空集合與是否有更改*/
	bool sym_check_naive() {
		bool terminal = true, nosym;
		pset p, q = set_copy(set_new(NUMINPUTS * 2), S);
		int i;

		for (i = 0, p = sym_table->data; i < NUMINPUTS; i++, p += sym_table->wsize) {
			nosym = setp_equal(p, cube.var_mask[i]);
			if (nosym) {
				var_remove(S, i);
			}
			else {
				terminal = false;
			}
		}

		state_is_change = false;
		if (!setp_equal(q, S)) {
			state_is_change = true;
		}
		FREE(q);

		return terminal;
	}

	/*回傳對稱配對表*/
	void return_Sym(pset_family A) {
		sf_copy(A, sym_table);
	}

	~Sym() {
		sf_free(sym_table);
		FREE(S);
	}
};

void sym_naive(pset_family F, pset_family R) {
	Sym sym(F, R);
	sym.main_naive();
}