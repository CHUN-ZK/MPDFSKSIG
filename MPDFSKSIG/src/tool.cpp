#include "header.h"

/*����pset��position����*/
char getvar(pset a, int position) {
	/*00 -> '?' 10 -> '0' 01 -> '1' 11 -> '2'*/
	return "?012"[GETINPUT(a, position)];
}

/*��Xpset�A�Ppbv1�A���L���׭���*/
string pstring(pset s, int n) {
	string s1(n, '1');

	for (int i = 0; i < n; i++) {
		s1[i] = is_in_set(s, i) ? '1' : '0';
	}
	return s1;
}

/*��Xpset_family�A���ΤFpstring*/
void sf_bm_printf(pset_family A) {
	pset p;
	int i;
	foreachi_set(A, i, p)
		cout << pstring(p, NUMINPUTS * 2) << endl;
	//printf("%s\n", i, pstring(p, A->sf_size).c_str());
	printf("\n");
}

/*�ˬd���pset��output�O�_�ۥ�*/
bool output_intersect(pset a, pset b) {
	for (int i = fow; i <= low; i++) {
		if (a[i] & b[i] & cube.mv_mask[i])
			return true;
	}
	return false;
}

/*���Ω�support-set�b���}��X�W*/
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

/*�]�wcube���������*/
void set_cube(int in, int out) {
	//save_cube_struct();
	cube.num_binary_vars = in;
	cube.num_vars = in + 1;
	cube.part_size = ALLOC(int, cube.num_vars);
	cube.part_size[cube.num_vars - 1] = out;
	cube_setup();
}

/*�Npset_family�����ǥ���*/
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

/*�^��pset��Don't care��mask*/
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

/*�^��pset���DDon't care��mask*/
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

/*�p��pset��Don't care���ƶq*/
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

/*�p��pset��Don't care���ƶq�A�æ^��Don't care��mask*/
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

/*�p����pset���ܼƦ�Don't care���ƶq*/
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

/*�p����pset���ܼƳ���Don't care���ƶq*/
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

/*�զX�ƾ�*/
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

/*��XTABLE�榡*/
void print_TABLE(TABLE t) {
	for (int i = 0; i < t.size(); i++) {
		for (int j = 0; j < t[i].distx.size(); j++) {
			//cout << t[i].distx[j].get_value() << " ";
			cout << t[i].distx[j] << " ";
		}
		cout << endl;
	}
}

/*����ƱƧ�pset_family*/
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

/*���P���󪺶Z���ˬd*/
/*�p����pset�S���ۥ�^��1�A�Ϥ�0*/
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

/*�p����pset�Z����1�^��1�A�Ϥ�0*/
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

/*�p����pset���Z���A�j��2�^��3�A�p�󵥩�2�|�b�^�ǶZ�����ͪ��ܼƦ�m*/
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

/*�p����pset���Z���A�j��2�^��3*/
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

/*�p����pset���Z��*/
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

/*�p���ƪ�support set�æ^��*/
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