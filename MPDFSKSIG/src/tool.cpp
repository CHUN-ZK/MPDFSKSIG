#include "header.h"

/*拿取pset中position的值*/
char getvar(pset a, int position) {
	/*00 -> '?' 10 -> '0' 01 -> '1' 11 -> '2'*/
	return "?012"[GETINPUT(a, position)];
}

/*輸出pset，同pbv1，但無長度限制*/
string pstring(pset s, int n) {
	string s1(n, '1');

	for (int i = 0; i < n; i++) {
		s1[i] = is_in_set(s, i) ? '1' : '0';
	}
	return s1;
}

/*輸出pset_family，應用了pstring*/
void sf_bm_printf(pset_family A) {
	pset p;
	int i;
	foreachi_set(A, i, p)
		cout << pstring(p, NUMINPUTS * 2) << endl;
	//printf("%s\n", i, pstring(p, A->sf_size).c_str());
	printf("\n");
}

/*檢查兩個pset的output是否相交*/
bool output_intersect(pset a, pset b) {
	for (int i = fow; i <= low; i++) {
		if (a[i] & b[i] & cube.mv_mask[i])
			return true;
	}
	return false;
}

/*應用於support-set在分開輸出上*/
pcover sep_sup_output(pset_family T, int i, pset sup) {
	pcover T1;
	pcube p, last, pdest, mask;

	mask = cube.var_mask[cube.output];
	T1 = new_cover(T->count);
	foreach_set(T, last, p) {
		if (is_in_set(p, i)) {
			pdest = GETSET(T1, T1->count++);
			INLINEset_and(pdest, p, sup);
			INLINEset_or(pdest, pdest, mask);
			RESET(pdest, PRIME);
		}
	}
	return T1;
}

/*設定cube的內部資料*/
void set_cube(int in, int out) {
	//save_cube_struct();
	cube.num_binary_vars = in;
	cube.num_vars = in + 1;
	cube.part_size = ALLOC(int, cube.num_vars);
	cube.part_size[cube.num_vars - 1] = out;
	cube_setup();
}

/*將pset_family的順序打亂*/
pset_family sf_rand(pset_family A) {
	pset_family p = sf_new(A->count, cube.size);
	int x;

	while (A->count > 0) {
		x = rand() % (A->count);
		sf_addset(p, GETSET(A, x));
		sf_delset(A, x);
	}
	sf_free(A);
	return p;
}

/*回傳pset中Don't care的mask*/
pset set_dc(pset r, pset a) {
	unsigned int x;
	int tmp;

	if ((tmp = cube.inword) != -1) {
		x = (a[tmp] >> 1) & a[tmp];
		x = x & cube.inmask;
		r[tmp] = x | (x << 1);
	}
	for (int i = 1; i < tmp; i++) {
		x = (a[i] >> 1) & a[i];
		x = x & DISJOINT;
		r[i] = x | (x << 1);
	}
	return r;
}

/*回傳pset中非Don't care的mask*/
pset set_ndc(pset r, pset a) {
	unsigned int x;
	int tmp;

	if ((tmp = cube.inword) != -1) {
		x = (a[tmp] >> 1) ^ a[tmp];
		x = x & cube.inmask;
		r[tmp] = x | (x << 1);
	}
	for (int i = 1; i < tmp; i++) {
		x = (a[i] >> 1) ^ a[i];
		x = x & DISJOINT;
		r[i] = x | (x << 1);
	}
	return r;
}

/*計算pset中Don't care的數量*/
int dc_count(pset a) {
	unsigned int x;
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		x = (a[tmp] >> 1) & a[tmp];
		x = x & cube.inmask;
		count += count_ones(x);
	}
	for (int i = 1; i < tmp; i++) {
		x = (a[i] >> 1) & a[i];
		x = x & DISJOINT;
		count += count_ones(x);
	}
	return count;
}

/*計算pset中Don't care的數量，並回傳Don't care的mask*/
int dc_count(pset r, pset a) {
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		r[tmp] = (a[tmp] >> 1) & a[tmp] & cube.inmask;
		count += count_ones(r[tmp]);
	}
	for (int i = 1; i < tmp; i++) {
		r[i] = (a[i] >> 1) & a[i] & DISJOINT;
		count += count_ones(r[i]);
	}
	return count;
}

/*計算兩個pset中變數有Don't care的數量*/
int set_ordc(pset a, pset b) {
	unsigned int x;
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		x = ((a[tmp] >> 1) & a[tmp]) | ((b[tmp] >> 1) & b[tmp]);
		x = x & cube.inmask;
		count += count_ones(x);
	}
	for (int i = 1; i < tmp; i++) {
		x = ((a[i] >> 1) & a[i]) | ((b[i] >> 1) & b[i]);
		x = x & DISJOINT;
		count += count_ones(x);
	}
	return count;
}

/*計算兩個pset中變數都為Don't care的數量*/
int set_anddc(pset a, pset b) {
	unsigned int x;
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		x = a[tmp] & b[tmp];
		x = (x >> 1) & x;
		x = x & cube.inmask;
		count += count_ones(x);
	}
	for (int i = 1; i < tmp; i++) {
		x = a[i] & b[i];
		x = (x >> 1) & x;
		x = x & DISJOINT;
		count += count_ones(x);
	}
	return count;
}

/*組合數學*/
cpp_int Cxy(int x, int y) {
	if (y < 0) return 0;
	else if (y == 0) return 1;

	int yy = y > (x - y) ? (x - y) : y;

	cpp_int xcount = 1, ycount = 1;
	for (int i = x; i > (x - yy); i--) {
		xcount *= i;
	}
	for (int i = yy; i > 0; i--) {
		ycount *= i;
	}

	return (xcount / ycount);
}

/*輸出TABLE格式*/
void print_TABLE(TABLE t) {
	for (int i = 0; i < t.size(); i++) {
		for (int j = 0; j < t[i].distx.size(); j++) {
			//cout << t[i].distx[j].get_value() << " ";
			cout << t[i].distx[j] << " ";
		}
		cout << endl;
	}
}

/*平行化排序pset_family*/
pset_family sf_parallel_sort(pset_family pf) {
	vector<pset> v(pf->count);
	int i;
	pset p;
	foreachi_set(pf, i, p) {
		pset tmp = set_new(pf->sf_size);
		tmp = set_copy(tmp, p);
		v[i] = tmp;
	}
	tbb::parallel_sort(v.begin(), v.end(), [&](const pset a, const pset b) {
		return set_ord(a) < set_ord(b);
		});

	pset_family R = sf_new(pf->count, pf->sf_size);
	R->count = pf->count;
	for (i = 0, p = R->data; i < v.size(); i++, p += R->wsize) {
		set_copy(p, v[i]);
	}
	sf_free(pf);

	return R;
}

/*不同條件的距離檢查*/
/*計算兩個pset沒有相交回傳1，反之0*/
bool notcdist(pset a, pset b) {
	unsigned int x;
	int tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j]; //10 ^ 01 = 11 => distance on this var
		x = (x >> 1) & x & DISJOINT;
		if (count_ones(x) > 0) {
			return false;
		}
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		if (count_ones(x) > 0) {
			return false;
		}
	}

	return true;
}

/*計算兩個pset距離為1回傳1，反之0*/
int cdist1(pset a, pset b) {
	unsigned int x;
	int tmp = cube.inword, ones; // save how many ones
	int dist = 0, n = 0x1, var = -1;
	bool dist_is_one = true;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		ones = count_ones(x);

		if (dist + ones > 1) return -1;
		else {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					var = (j - 1) * varnum + i / 2;
					ones--;
				}
			}
		}
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		ones = count_ones(x);
		if (dist + ones > 1) return -1;
		else if (ones > 0) {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					var = (tmp - 1) * varnum + i / 2;
					ones--;
				}
			}
		}
	}
	return var;
}

/*計算兩個pset的距離，大於2回傳3，小於等於2會在回傳距離產生的變數位置*/
int cdist123(pset a, pset b, int dist_on[2]) {
	unsigned int x;
	int tmp, ones; // save how many ones
	int counts = 0; // save dist's position
	int dist = 0, n = 0x1;

	tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		ones = count_ones(x);
		if (dist + ones > 2) return 3;
		else if (ones > 0) {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					dist_on[counts] = (j - 1) * varnum + i / 2;
					ones--;
					counts++;
				}
			}
		}
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		ones = count_ones(x);
		if (dist + ones > 2) return 3;
		else if (ones > 0) {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					dist_on[counts] = (tmp - 1) * varnum + i / 2;
					ones--;
					counts++;
				}
			}
		}
	}

	return dist;
}

/*計算兩個pset的距離，大於2回傳3*/
int cdist123(pset a, pset b) {
	unsigned int x;
	int tmp; // save how many ones
	int dist = 0;

	tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		dist += count_ones(x);
		if (dist > 2) return 3;
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		dist += count_ones(x);
		if (dist > 2) return 3;
	}

	return dist;
}

/*計算兩個pset的距離*/
int cdistN(pset a, pset b) {
	unsigned int x;
	int tmp;
	int dist = 0;

	tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		dist += count_ones(x);
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		dist += count_ones(x);
	}

	return dist;
}

/*計算函數的support set並回傳*/
pset_family support_set(pset_family F, pset_family R) {
	pset_family supp = sf_new(NUMOUTPUTS, cube.size);
	int fnum, rnum, var, a, b;

	unsigned int x;
	pset fp, rp, k;

	for (int i = 0; i < NUMOUTPUTS; i++) {
		sf_addset(supp, set_new(cube.size));
		k = supp->data + supp->wsize * (i);
		set_insert(k, NUMINPUTS * 2 + i);
	}

	a = BPI - NUMINPUTS * 2 % BPI;
	b = NUMOUTPUTS - (a + (low - fow - 1) * BPI);

	foreachi_set(F, fnum, fp) {
		foreachi_set(R, rnum, rp) {
			if (output_intersect(fp, rp)) {
				var = cdist1(fp, rp);
				if (var == -1) continue;

				x = fp[fow] & rp[fow] & cube.mv_mask[fow];
				if (x) {
					for (int i = BPI - a, n = 0x1 << i; i < BPI; i++, n <<= 1) {
						if (x & n) {
							k = supp->data + supp->wsize * (i - (BPI - a));
							var_insert(k, var, 2);
						}
					}
				}

				if (low != fow) {
					for (int i = fow + 1; i < low; i++) {
						x = fp[i] & rp[i] & cube.mv_mask[i];
						if (x) {
							for (int j = 0, n = 0x1 << j; j < BPI; j++, n <<= 1) {
								if (x & n) {
									k = supp->data + supp->wsize * (a + (i - fow - 1) * BPI + j);
									var_insert(k, var, 2);
								}
							}
						}
					}
					x = fp[low] & rp[low] & cube.mv_mask[low];
					if (x) {
						for (int i = 0, n = 0x1 << i; i < b; i++, n <<= 1) {
							if (x & n) {
								k = supp->data + supp->wsize * (a + (low - fow - 1) * BPI + i);
								var_insert(k, var, 2);
							}
						}
					}
				}
			}
		}
	}
	return supp;
}