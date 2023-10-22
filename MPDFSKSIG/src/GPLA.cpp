#include <random>
#include "header.h"

/*
	生成大型單輸出電路
*/

/*設定cube大小與Random種子*/
/*cube.size = input_size * 2 + 1*/
void general_setup(int input, int output) {
	set_cube(input, output);
	srand(time(NULL));
}

/*min <= randdata(max, min) <= max*/
int randdata(int max, int min) {
	return rand() % (max - min + 1) + min;
}

/*input_size * 0.2  <  number of don't care pset  < input_size * 0.8 */
pset randpset(pset p, int size, int dc) {
	p = set_clear(p, size * 2 + 1);
	set_insert(p, size * 2);
	
	int dcpart = dc;
	int i, x, pos;
	for (i = 0; i < dcpart; i++) {
		pos = randdata(size - 1, 0);
		var_insert(p, pos, 2);
	}
	for (i = 0; i < size; i++) {
		if(!is_in_set(p, i * 2))
			var_insert(p, i, rand() % 2);
	}
	return p;
}

/*隨機生成don't care數量為high_p至low_p之間的pset*/
pset randpset_p(pset p, int size, double high_p, double low_p) {
	p = set_clear(p, size * 2 + 1);
	set_insert(p, size * 2);

	int dcpart = randdata(high_p * size, low_p * size);
	int i, x, pos;
	for (i = 0; i < dcpart; i++) {
		pos = randdata(size - 1, 0);
		var_insert(p, pos, 2);
	}
	for (i = 0; i < size; i++) {
		var_insert(p, i, rand() % 2);
	}
	return p;
}

/*檢查A中所有積項是否包含p，false為有包含*/
bool check_imply(pset_family A, pset p) {
	pset k;
	int i;

	foreachi_set(A, i, k) {
		/*if distance of k and p = 0, non feasible*/
		if (notcdist(k, p))
			return false;
	}
	return true;
}

/*全隨機生成大型電路*/
void GeneralPLA_allrandom(pset_family& F, pset_family& R, int input_size, int onset_size, int offset_size, int on_off) {

	pset p = set_new(input_size * 2 + 1), k;
	int while_time = 100, is_time_to_end = 0;
	int i, pos, pos2, dc;
	char c;

	/*step 1 onset is generated to the specified number, must control the amount of don't care*/
	for (i = F->count; i < onset_size; i++) {
		dc = 10; //input_size / 3
		p = randpset(p, input_size, dc);

		/*check distance of new pset and R */
		if (check_imply(R, p)) {
			sf_addset(F, p);
			is_time_to_end = 0;
		}
		else {
			i--;
		}
		is_time_to_end++;
		if (is_time_to_end == while_time) {
			break;
		}
	}

	/*step 2 Generate a specified number of offsets, requiring a distance of 1 or 2 from the onset(on-off), and record the removed pair*/
	is_time_to_end = 0;
	for (i = R->count; i < on_off; i++) {
		/*choose which pset of F*/
		pos = randdata(F->count - 1, 0);
		k = F->data + F->wsize * pos;
		p = set_copy(p, k);
		/*choose which variable of pset, this variable can't be don't care*/
		if (rand() % 3 > 0) {
			pos = randdata(input_size - 1, 0);
			c = getvar(k, pos);
			while (c == '2') {
				pos = randdata(input_size - 1, 0);
				c = getvar(k, pos);
			}
			if (c == '0') {
				set_remove(p, pos * 2);
				set_insert(p, pos * 2 + 1);
			}
			else {
				set_remove(p, pos * 2 + 1);
				set_insert(p, pos * 2);
			}
		}
		else {
			pos = randdata(input_size - 1, 0);
			c = getvar(k, pos);
			while (c == '2') {
				pos = randdata(input_size - 1, 0);
				c = getvar(k, pos);
			}
			if (c == '0') {
				set_remove(p, pos * 2);
				set_insert(p, pos * 2 + 1);
			}
			else {
				set_remove(p, pos * 2 + 1);
				set_insert(p, pos * 2);
			}
			pos2 = randdata(input_size - 1, 0);
			c = getvar(k, pos2);
			while (c == '2' && pos != pos2) {
				pos2 = randdata(input_size - 1, 0);
				c = getvar(k, pos);
			}
			if (c == '0') {
				set_remove(p, pos2 * 2);
				set_insert(p, pos2 * 2 + 1);
			}
			else {
				set_remove(p, pos2 * 2 + 1);
				set_insert(p, pos2 * 2);
			}
		}

		/*check distance of new pset and F */
		if (check_imply(F, p)) {
			sf_addset(R, p);
			is_time_to_end = 0;
		}
		else {
			i--;
		}
		is_time_to_end++;
		if (is_time_to_end == while_time) {
			break;
		}
	}

	/*step 3 onset generated to the specified number*/ // redundant if doesn't using implies 
	is_time_to_end = 0;
	for (i = F->count; i < onset_size; i++) {
		dc = 10; //input_size / 3
		p = randpset(p, input_size, dc);
		/*check distance of new pset and R */
		if (check_imply(R, p)) {
			sf_addset(F, p);
			is_time_to_end = 0;
		}
		else {
			i--;
		}
		is_time_to_end++;
		if (is_time_to_end == while_time) {
			break;
		}
	}

	/*step 4-1 offset generated to the specified number*/		
	is_time_to_end = 0;
	for (i = R->count; i < offset_size; i++) {
		p = randpset(p, input_size, 10);
		/*check distance of new pset and F */
		if (check_imply(F, p)) {
			sf_addset(R, p);
			is_time_to_end = 0;
		}
		else {
			i--;
		}
		is_time_to_end++;
		if (is_time_to_end == while_time) {
			break;
		}
	}

	/*step 下下策，使用已經在pset_family中的成員進行填充*/
	for (i = F->count; i < onset_size; i++) {
		pos = randdata(F->count-1, 0);
		k = F->data + F->wsize * pos;
		sf_addset(F, k);
	}
	for (i = R->count; i < offset_size; i++) {
		pos = randdata(R->count-1, 0);
		k = R->data + R->wsize * pos;
		sf_addset(R, k);
	}
	set_free(p);
}

/*生成沒有函數對稱性電路*/
bool GeneralPLA_nonSym(pset_family& F, pset_family& R, int input_size, int onset_size, int offset_size, int on_off) {
	pset p = set_new(input_size * 2 + 1);
	pset q = set_new(input_size * 2 + 1);
	pset p1 = set_new(input_size * 2 + 1);
	pset q1 = set_new(input_size * 2 + 1);
	int var, rvc, x, dc, i, j; //remove variable count
	ullg while_time = 1000000000, is_time_to_end = 0;
	boolean e = false;
	bool ok = false;
	bool fok, rok, fok1, rok1;

	/*remove symmetry table of var*/
	for (var = 0; var < input_size - 1 && !e; var++) {
		ok = false;
		rvc = var + 1; //remove from var + 1

		while (!ok && !e) {
			if (rand() % 3 == 0) {
				bool ne_remove = 0;
				/*NE & E*/
				is_time_to_end = 0;
				do {
					p = set_clear(p, input_size * 2 + 1);
					set_insert(p, input_size * 2);
					q = set_clear(q, input_size * 2 + 1);
					set_insert(q, input_size * 2);

					p1 = set_clear(p1, input_size * 2 + 1);
					set_insert(p1, input_size * 2);
					q1 = set_clear(q1, input_size * 2 + 1);
					set_insert(q1, input_size * 2);

					x = rand() % 2; //0 or 1
					i = x;
					var_insert(p, var, x);
					var_insert(p, rvc, (1 - x));
					var_insert(q, var, (1 - x));
					var_insert(q, rvc, x);
					var_insert(p1, var, x);
					var_insert(p1, rvc, x);
					var_insert(q1, var, (1 - x));
					var_insert(q1, rvc, (1 - x));
					for (j = 0; j < input_size; j++) {
						if (j != var && j != rvc) {
							x = rand() % 3;
							var_insert(p, j, x);
							var_insert(q, j, x);
							var_insert(p1, j, x);
							var_insert(q1, j, x);
						}
					}
					fok = check_imply(R, p);
					fok1 = check_imply(R, p1);
					rok = check_imply(F, q);
					rok1 = check_imply(F, q1);
					is_time_to_end++;

				} while (!fok && !rok && !fok1 && !rok1 && is_time_to_end < while_time);
				if (fok && rok && fok1 && rok1) {
					sf_addset(F, p);
					sf_addset(F, p1);
					sf_addset(R, q);
					sf_addset(R, q1);
				}
				else {
					e = true;
				}
			}
			else {
				/*don't make too much don't care*/
				dc = randdata(input_size/10, 1);
				if (dc + rvc > input_size) {
					dc = input_size - rvc;
				}
				is_time_to_end = 0;
				do {
					p = set_clear(p, input_size * 2 + 1);
					set_insert(p, input_size * 2);
					q = set_clear(q, input_size * 2 + 1);
					set_insert(q, input_size * 2);


					for (j = 0; j < input_size; j++) {
						if (j == var) {
							x = rand() % 2;
							var_insert(p, j, x);
							var_insert(q, j, (1 - x));
						}
						else if (j >= rvc && j < rvc + dc) {
							var_insert(p, j, 2);
							var_insert(q, j, 2);
						}
						else {
							x = rand() % 3;
							var_insert(p, j, x);
							var_insert(q, j, x);
						}
					}
					fok = check_imply(F, q);
					rok = check_imply(R, p);
					is_time_to_end++;

				} while (!fok && !rok && is_time_to_end < while_time);
				if (fok && rok) {
					sf_addset(F, p);
					sf_addset(R, q);
					rvc += dc;
					if (rvc == input_size) {
						ok = true;
					}
				}
				else {
					e = true;
				}
			}
		}


	}

	set_free(p);
	set_free(q);
	set_free(p1);
	set_free(q1);
	if (e) {
		return false;
	}
	return true;
}
